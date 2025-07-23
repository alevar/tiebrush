use crate::samreader::SAMReader;
use crate::commons::*;

use std::fs::File;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use clap::{Args, ArgGroup};
use std::collections::{HashMap, hash_map::Entry};
use rust_htslib::bam::{Record,record::{Aux}, Header, HeaderView};
use rust_htslib::bam::ext::BamRecordExtensions;

#[derive(Args, Debug)]
#[command(group(
    ArgGroup::new("inputs")
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
}

/// JuncMat is a matrix of junction data for a single (seqid, strand).
/// user must assert that the seqid and strand are the same for all entries in the matrix.
/// otherwise errors will be thrown.
#[derive(Debug)]
pub struct JuncMat {
    pub seqid: i32,
    data: HashMap<(Strand,i64, i64), u64>,
}

impl JuncMat {
    pub fn new(seqid: i32) -> Self {
        Self { seqid, data: HashMap::new() }
    }

    pub fn insert(&mut self, strand: Strand, start: i64, end: i64, val: u64) {
        self.data.insert((strand, start, end), val);
    }

    pub fn get(&self, strand: Strand, start: i64, end: i64) -> Option<&u64> {
        self.data.get(&(strand, start, end))
    }

    pub fn get_mut(&mut self, strand: Strand, start: i64, end: i64) -> Option<&mut u64> {
        self.data.get_mut(&(strand, start, end))
    }

    /// Increment the value for a given junction, or insert with 1 if not present (for counts).
    pub fn increment(&mut self, strand: Strand, start: i64, end: i64)
    {
        match self.data.entry((strand, start, end)) {
            Entry::Occupied(mut e) => *e.get_mut() += 1,
            Entry::Vacant(e) => { e.insert(1); },
        }
    }
}

impl Default for JuncMat {
    fn default() -> Self {
        Self {
            seqid: 0,
            data: HashMap::new(),
        }
    }
}

#[derive(Debug)]
pub struct TBCov {
    seqid: i32,
    start: i64,
    end: i64,
    cov: Vec<i32>,
    junc: JuncMat,
    
    store_cov: bool,
    store_junc: bool,
}

impl TBCov {
    pub fn new(seqid: i32, start: i64, end: i64, val: i32) -> Self {
        Self {  seqid,
                start, end,
                cov: vec![val as i32; (end - start) as usize],
                junc: JuncMat::default(),
                store_cov: true,
                store_junc: true,
        }
    }

    pub fn new_from_record(record: &Record) -> Self {
        let start = record.pos();
        let end = record.reference_end();
        Self::new(record.tid(), start, end, 1)
    }

    pub fn set_store_cov(&mut self, store_cov: bool) {
        self.store_cov = store_cov;
    }

    pub fn set_store_junc(&mut self, store_junc: bool) {
        self.store_junc = store_junc;
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
    pub fn add_record(&mut self, record: &Record) -> anyhow::Result<()> {
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
            None => 1,
        };
        if yc > 1 {
            println!("yc: {}", yc);
        }

        if record.reference_end() > self.end {
            // extend the cov vec
            let extend_len = (record.reference_end() - self.end) as usize;
            self.cov.extend(vec![0; extend_len]);
            self.end = record.reference_end();
        }

        for pos in record.reference_positions() {
            let vec_pos = (pos - self.start) as usize;
            self.cov[vec_pos] += yc;
        }
        
        Ok(())
    }
}

impl Default for TBCov {
    fn default() -> Self {
        Self::new(0, 0, 0, 0)
    }
}

pub struct CovCMD {
    cov_args: CovArgs,
    cov_writer: Option<BufWriter<File>>,
    junc_writer: Option<BufWriter<File>>,
}

impl CovCMD {
    pub fn new(args: CovArgs) -> anyhow::Result<Self> {
        // verify all provided alignments exists
        for alignment in &args.input_alignments {
            if !alignment.exists() {
                anyhow::bail!("Alignment file {} does not exist", alignment.display());
            }
        }

        let cov_writer = match args.coverage {
            Some(ref cov_path) => {
                let mut writer = BufWriter::new(File::create(cov_path)?);
                writeln!(writer, "track type=bedGraph")?;
                Some(writer)
            },
            None => None,
        };
        let junc_writer = match args.junctions {
            Some(ref junc_path) => {
                Some(BufWriter::new(File::create(junc_path)?))
            },
            None => None,
        };

        Ok(Self {
            cov_args: args,
            cov_writer,
            junc_writer,
        })
    }

    pub fn flush_cov(&mut self, tbcov: &TBCov, header: &HeaderView) -> anyhow::Result<()> {
        match self.cov_writer {
            Some(ref mut writer) => {
                if tbcov.cov.is_empty() {
                    return Ok(());
                }
                
                let mut first_pos = tbcov.start;
                let mut last_pos = tbcov.start;
                let mut pos_cov = tbcov.cov[0];
                for (i, val) in tbcov.cov.iter().enumerate() {
                    if *val == 0 {
                        continue;
                    }
                    let pos = tbcov.start + i as i64;
                    if *val != pos_cov {
                        let seqname = std::str::from_utf8(header.tid2name(tbcov.seqid as u32)).unwrap();
                        writeln!(writer, "{}\t{}\t{}\t{}", seqname, first_pos, last_pos + 1, pos_cov)?;
                        first_pos = pos;
                        pos_cov = *val;
                    }
                    last_pos = pos;
                }
                // Write the last interval
                let seqname = std::str::from_utf8(header.tid2name(tbcov.seqid as u32)).unwrap();
                writeln!(writer, "{}\t{}\t{}\t{}", seqname, first_pos, last_pos + 1, pos_cov)?;
            }
            None => {}
        }
        Ok(())
    }

    pub fn flush_junc(&self, tbcov: &TBCov, header: &HeaderView) -> anyhow::Result<()> {
        match self.cov_args.junctions {
            Some(ref junc_path) => {
                // println!("junc: {:?}", tbcov.junc);
            }
            None => {}
        }
        Ok(())
    }

    pub fn run(&mut self) -> anyhow::Result<()> {
        let mut sam_reader = SAMReader::new(&self.cov_args.input_alignments)?;
        let header = sam_reader.get_header();
        let header_view = HeaderView::from_header(header);
        // get first record;
        let mut tbcov = match sam_reader.next() {
            Some(record) => TBCov::new_from_record(&record),
            None => {
                anyhow::bail!("No records found in the input alignments");
            }
        };

        // need to somehow map from seqid to seqname 
        
        for record in sam_reader {
            if record.is_unmapped() {
                continue;
            }

            if record.tid() != tbcov.seqid || record.pos() > tbcov.end {
                self.flush_cov(&tbcov, &header_view)?;
                self.flush_junc(&tbcov, &header_view)?;
                tbcov = TBCov::new_from_record(&record);
            }
            else {
                tbcov.add_record(&record)?;
            }
        }
        self.flush_cov(&tbcov, &header_view)?;
        self.flush_junc(&tbcov, &header_view)?;

        Ok(())
    }
}

pub fn run(args: CovArgs) -> anyhow::Result<()> {
    let mut cmd = CovCMD::new(args)?;
    cmd.run()?;

    Ok(())
}