use crate::samreader::TBSAMReader;
use crate::commons::*;

use std::fs::File;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use clap::{Args, ArgGroup};
use std::collections::{HashMap, hash_map::Entry};
use rust_htslib::bam::{Record, HeaderView, record::Cigar, ext::BamRecordExtensions};
use rust_htslib::htslib;

use bigtools::{BigWigWrite, Value};
use bigtools::beddata::BedParserStreamingIterator;

use tokio::runtime::Runtime;

#[derive(Args, Debug)]
#[command(group(
    ArgGroup::new("outputs")
    .required(true)
    .multiple(true)
        .args(["coverage", "junctions"])
))]
pub struct CovArgs {
    /// BedGraph (or BedWig with '-W') file with coverage for all mapped bases.
    #[arg(short, long)]
    pub coverage: Option<PathBuf>,
    /// BED file with coverage of all splice-junctions in the input file.
    #[arg(short, long)]
    pub junctions: Option<PathBuf>,
    /// Save coverage in BigWig format. Default output is in Bed format.
    #[arg(short='W', long)]
    pub bigwig: bool,
    /// One or more alignment files in SAM/BAM/CRAM format.
    #[arg(required = true)]
    pub input_alignments: Vec<PathBuf>,
    /// Exclude supplementary alignments. Supplementary alignments are included by default.
    #[arg(short='S', long)]
    pub exclude_supplementary: bool,
    /// Include only flags. Reads with the provided bits set in the flag field are included.
    #[arg(short='f', long)]
    pub include_flags: Option<u16>,
    /// Exclude flags. Reads with the provided bits set in the flag field are excluded.
    #[arg(short='F', long)]
    pub exclude_flags: Option<u16>,
}

/// JuncMat is a matrix of junction data for a single (seqid, strand).
/// user must assert that the seqid and strand are the same for all entries in the matrix.
/// otherwise errors will be thrown.
#[derive(Debug)]
pub struct JuncMat {
    data: HashMap<(Strand,i64, i64), f64>,
}

impl JuncMat {
    pub fn new() -> Self {
        Self {data: HashMap::new() }
    }

    pub fn insert(&mut self, strand: Strand, start: i64, end: i64, val: f64) {
        self.data.insert((strand, start, end), val);
    }

    pub fn get(&self, strand: Strand, start: i64, end: i64) -> Option<&f64> {
        self.data.get(&(strand, start, end))
    }

    pub fn get_mut(&mut self, strand: Strand, start: i64, end: i64) -> Option<&mut f64> {
        self.data.get_mut(&(strand, start, end))
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Increment the value for a given junction, or insert with the provided value if not present (for counts).
    pub fn increment(&mut self, strand: Strand, start: i64, end: i64, val: f64)
    {
        match self.data.entry((strand, start, end)) {
            Entry::Occupied(mut e) => *e.get_mut() += val,
            Entry::Vacant(e) => { e.insert(val); },
        }
    }
}

impl Default for JuncMat {
    fn default() -> Self {
        Self {
            data: HashMap::new(),
        }
    }
}

#[derive(Debug)]
pub struct TBCov {
    seqid: i32,
    start: i64,
    end: i64,
    cov: Vec<f64>,
    junc: JuncMat,
}

impl TBCov {
    pub fn new(seqid: i32, start: i64, end: i64, val: f64) -> Self {
        Self {  seqid,
                start, end,
                cov: vec![val as f64; (end - start) as usize],
                junc: JuncMat::default(),
        }
    }

    pub fn new_from_record(record: &Record, store_cov: bool, store_junc: bool) -> Self {
        let start = record.pos();
        let end = record.reference_end();
        let mut tbcov = Self::new(record.tid(), start, end, 0.0);
        tbcov.add_record(record, store_cov, store_junc);
        tbcov
    }
    
    /// Adds a record's coverage to the BCov region.
    ///
    /// # Errors
    /// - Returns an error if the record's `tid` (reference sequence ID) does not match this BCov's `seqid`.
    /// - Returns an error if the record start preceeds the start of the BCov.
    ///
    /// # Note
    /// - The caller is responsible for ensuring only compatible records are passed to this method.
    /// - Unmapped records are ignored (no error is returned).
    pub fn add_record(&mut self, record: &Record, store_cov: bool, store_junc: bool) -> anyhow::Result<()> {
        if record.tid() != self.seqid {
            anyhow::bail!("Record tid {} does not match BCov seqid {}", record.tid(), self.seqid);
        }
        if record.is_unmapped() {
            return Ok(());
        }
        if record.pos() < self.start {
            anyhow::bail!("Record pos {} is before BCov start {}. Reads must be sorted by coordinate.", record.pos(), self.start);
        }

        // get YC tag if available as coverage otherwise use 1
        let yc = match get_yc_tag(&record)? {
            Some(yc) => yc,
            None => 1.0,
        };

        if record.reference_end() > self.end && store_cov {
            // extend the cov vec
            let extend_len = (record.reference_end() - self.end) as usize;
            self.cov.extend(vec![0.0; extend_len]);
            self.end = record.reference_end();
        }

        // iterate using cigar operations
        let mut ref_pos = record.pos();
        for op in record.cigar().iter() {
            match op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    if store_cov {
                        for _ in 0..*len {
                            self.cov[(ref_pos - self.start) as usize] += yc;
                            ref_pos += 1;
                        }
                    }
                    else {
                        ref_pos += *len as i64;
                        continue;
                    }
                }
                Cigar::Ins(_len) | Cigar::SoftClip(_len) => {
                }
                Cigar::Del(len) => {
                    ref_pos += *len as i64;
                }
                Cigar::RefSkip(len) => {
                    // add to junc mat
                    if store_junc {
                        let strand = get_strand(&record)?;
                        self.junc.increment(strand, ref_pos, ref_pos + *len as i64, yc.try_into().unwrap());
                    }
                    ref_pos += *len as i64;
                }
                Cigar::HardClip(_) => {} // no advance
                Cigar::Pad(_) => panic!("Padding (Cigar::Pad) is not supported."), //padding is only used for multiple sequence alignment
            }
        }
        
        Ok(())
    }
}

impl Default for TBCov {
    fn default() -> Self {
        Self::new(0, 0, 0, 0.0)
    }
}

trait CoverageWriter {
    fn write_interval(&mut self, chrom: &str, start: i64, end: i64, value: f32) -> anyhow::Result<()>;
    fn flush(&mut self) -> anyhow::Result<()>;
}

struct BedGraphWriter {
    writer: BufWriter<File>,
}

impl BedGraphWriter {
    fn new(path: PathBuf) -> anyhow::Result<Self> {
        let mut writer = BufWriter::new(File::create(path)?);
        writeln!(writer, "track type=bedGraph")?;
        Ok(Self { writer })
    }
}

impl CoverageWriter for BedGraphWriter {
    fn write_interval(&mut self, chrom: &str, start: i64, end: i64, value: f32) -> anyhow::Result<()> {
        writeln!(self.writer, "{}\t{}\t{}\t{}", chrom, start, end, value)?;
        Ok(())
    }

    fn flush(&mut self) -> anyhow::Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

struct BigWigWriter {
    writer: Option<BigWigWrite<File>>,
    intervals: Vec<(String, i64, i64, f32)>,
}

impl BigWigWriter {
    fn new(mut path: PathBuf, header: &HeaderView) -> anyhow::Result<Self> {
        // check path and add .bw if not present
        let ext = path.extension().unwrap_or_default().to_string_lossy().to_lowercase();
        if ext != "bw" && ext != "bigwig" {
            path.set_extension("bw");
        }

        // Build chromosome map from BAM header
        let mut chrom_map = HashMap::new();
        for i in 0..header.target_count() {
            let name = std::str::from_utf8(header.tid2name(i))
                .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in chromosome name: {}", e))?;
            let length = header.target_len(i).unwrap_or(0) as u32;
            chrom_map.insert(name.to_string(), length);
        }

        let writer = BigWigWrite::create_file(path, chrom_map)?;
        Ok(Self {
            writer: Some(writer),
            intervals: Vec::new(),
        })
    }
}

impl CoverageWriter for BigWigWriter {
    fn write_interval(&mut self, chrom: &str, start: i64, end: i64, value: f32) -> anyhow::Result<()> {
        self.intervals.push((chrom.to_string(), start, end, value));
        Ok(())
    }

    fn flush(&mut self) -> anyhow::Result<()> {
        if !self.intervals.is_empty() {
            // Sort intervals by chromosome and position
            self.intervals.sort_by(|a, b| {
                a.0.cmp(&b.0).then(a.1.cmp(&b.1))
            });
            
            let data_intervals: Vec<_> = self.intervals.iter().map(|(chrom, start, end, value)| {
                (chrom.as_str(), Value {
                    start: *start as u32,
                    end: *end as u32,
                    value: *value,
                })
            }).collect();

            let data = BedParserStreamingIterator::wrap_infallible_iter(
                data_intervals.into_iter(),
                true // sorted
            );

            if let Some(writer) = self.writer.take() {
                let runtime = Runtime::new()?;
                let _ = writer.write(data, runtime);
            }
            
            self.intervals.clear();
        }
        Ok(())
    }
}

pub struct CovCMD {
    cov_args: CovArgs,
    cov_writer: Option<Box<dyn CoverageWriter>>,
    junc_writer: Option<BufWriter<File>>,
    junc_counter: u64, // used for naming junctions
    sam_reader: TBSAMReader,
    header_view: HeaderView,
}

impl CovCMD {
    pub fn new(args: CovArgs) -> anyhow::Result<Self> {
        // verify all provided alignments exists
        for alignment in &args.input_alignments {
            if !alignment.exists() {
                anyhow::bail!("Alignment file {} does not exist", alignment.display());
            }
        }

        // Create SAM reader to get header for BigWig initialization
        let req_fields = htslib::sam_fields_SAM_QNAME | htslib::sam_fields_SAM_FLAG | htslib::sam_fields_SAM_RNAME | htslib::sam_fields_SAM_POS | htslib::sam_fields_SAM_CIGAR | htslib::sam_fields_SAM_AUX;
        let sam_reader = TBSAMReader::new_with_fields(&args.input_alignments, req_fields)?;
        let header = sam_reader.get_header();
        let header_view = HeaderView::from_header(header);

        let cov_writer: Option<Box<dyn CoverageWriter>> = match args.coverage {
            Some(ref cov_path) => {
                if args.bigwig {
                    Some(Box::new(BigWigWriter::new(cov_path.clone(), &header_view)?))
                } else {
                    Some(Box::new(BedGraphWriter::new(cov_path.clone())?))
                }
            },
            None => None,
        };
        let junc_writer = match args.junctions {
            Some(ref junc_path) => {
                let mut writer = BufWriter::new(File::create(junc_path)?);
                writeln!(writer, "track name=junctions")?;
                Some(writer)
            },
            None => None,
        };

        Ok(Self {
            cov_args: args,
            cov_writer,
            junc_writer,
            junc_counter: 1,
            sam_reader,
            header_view,
        })
    }

    pub fn flush_cov(&mut self, tbcov: &TBCov) -> anyhow::Result<()> {
        if tbcov.cov.is_empty() {
            return Ok(());
        }
        
        let seqname = std::str::from_utf8(self.header_view.tid2name(tbcov.seqid as u32))
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence name: {}", e))?
            .to_string();

        let mut first_pos = tbcov.start;
        let mut last_pos = tbcov.start;
        let mut pos_cov = tbcov.cov[0];
        for (i, val) in tbcov.cov.iter().enumerate() {
            if *val == 0.0 {
                continue;
            }
            let pos = tbcov.start + i as i64;
            if *val != pos_cov {
                self.cov_writer.as_mut().unwrap().write_interval(&seqname, first_pos, last_pos + 1, pos_cov as f32)?;
                first_pos = pos;
                pos_cov = *val;
            }
            last_pos = pos;
        }
        // Write the last interval
        self.cov_writer.as_mut().unwrap().write_interval(&seqname, first_pos, last_pos + 1, pos_cov as f32)?;
        
        self.cov_writer.as_mut().unwrap().flush()?;
        Ok(())
    }

    pub fn flush_junc(&mut self, tbcov: &TBCov) -> anyhow::Result<()> {
        match self.junc_writer {
            Some(ref mut writer) => {
                if tbcov.junc.is_empty() {
                    return Ok(());
                }

                let mut entries: Vec<_> = tbcov.junc.data.iter().collect();
                entries.sort_by(|a, b| {
                    let ((_, start_a, end_a), _) = a;
                    let ((_, start_b, end_b), _) = b;
                    start_a.cmp(start_b).then(end_a.cmp(end_b))
                });

                for ((strand, start, end), count) in entries {
                    let seqname = std::str::from_utf8(self.header_view.tid2name(tbcov.seqid as u32)).unwrap();
                    let junc_name = format!("JUNC{:08}", self.junc_counter);
                    writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}", seqname, start, end, junc_name, count, strand)?;
                    self.junc_counter += 1;
                }
            }
            None => {}
        }
        Ok(())
    }

    pub fn run(&mut self) -> anyhow::Result<()> {
        // get first record;
        let mut tbcov = match self.sam_reader.next() {
            Some(tb_record) => {
                TBCov::new_from_record(&tb_record.record(),self.cov_args.coverage.is_some(),self.cov_args.junctions.is_some())
            },
            None => {
                anyhow::bail!("No records found in the input alignments");
            }
        };

        // Process remaining records
        while let Some(tb_record) = self.sam_reader.next() {
            let record = tb_record.record();
            if record.is_unmapped() {
                continue;
            }

            // check flags
            if self.cov_args.exclude_supplementary && record.is_supplementary() {
                continue;
            }
            if self.cov_args.include_flags.is_some() && !flags_set(&record, self.cov_args.include_flags.unwrap()) {
                continue;
            }
            if self.cov_args.exclude_flags.is_some() && flags_set(&record, self.cov_args.exclude_flags.unwrap()) {
                continue;
            }

            if record.tid() != tbcov.seqid || record.pos() > tbcov.end {
                if self.cov_args.coverage.is_some() {
                    self.flush_cov(&tbcov)?;
                }
                if self.cov_args.junctions.is_some() {
                    self.flush_junc(&tbcov)?;
                }
                tbcov = TBCov::new_from_record(&record,self.cov_args.coverage.is_some(),self.cov_args.junctions.is_some());
            }
            else {
                tbcov.add_record(&record,self.cov_args.coverage.is_some(),self.cov_args.junctions.is_some())?;
            }
        }
        if self.cov_args.coverage.is_some() {
            self.flush_cov(&tbcov)?;
        }
        if self.cov_args.junctions.is_some() {
            self.flush_junc(&tbcov)?;
        }

        Ok(())
    }
}

pub fn run(args: CovArgs) -> anyhow::Result<()> {
    let mut cmd = CovCMD::new(args)?;
    cmd.run()?;

    Ok(())
}