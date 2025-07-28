use crate::samreader::{TBSAMReader, TBSAMReaderRecord};
use crate::commons::*;

use std::{path::PathBuf, vec};
use std::collections::HashMap;
use clap::{Args, ArgGroup};
use rust_htslib::bam::{Record, record::Cigar, Format, ext::BamRecordExtensions, Writer};
use rust_htslib::htslib;

#[derive(Debug, Default)]
pub struct TBStats {
    pub total_input_reads: u64,
    pub filtered_reads: u64,
    pub processed_reads: u64,
    
    pub total_output_reads: u64,
}

impl TBStats {
    fn compression_ratio(&self) -> f64 {
        if self.total_output_reads == 0 {
            0.0
        } else {
            (self.processed_reads as f64) / (self.total_output_reads as f64)
        }
    }
    
    fn print_summary(&self) {
        println!("=== TieBrush Deduplication Summary ===");
        println!("Total input reads: {}", self.total_input_reads);
        println!("Filtered reads: {}", self.filtered_reads);
        println!("Input reads processed: {}", self.processed_reads);
        println!("Output reads written: {}", self.total_output_reads);
        println!("Compression ratio: {:.2}x", self.compression_ratio());
    }
}

#[derive(Args, Debug)]
#[command(group(
    ArgGroup::new("cmp_group")
        .multiple(false)
        .args(["full", "clip", "exon"])
))]
pub struct BrushArgs {
    /// Output TieBrush file. Format is configured from the file extension. 
    /// Only BAM and SAM formats are currently supported for output.
    #[arg(short, long, required = true)]
    pub output: Option<PathBuf>,
    /// One or more alignment files in SAM/BAM/CRAM format.
    #[arg(required = true)]
    pub input_alignments: Vec<PathBuf>,
    /// If enabled, only reads with the same CIGAR and MD strings will be grouped and collapsed.
    #[arg(short='L', long="full")]
    pub full: bool,
    /// If enabled, reads will be grouped by clipped CIGAR string. 
    /// In this mode 5S10M5S and 3S10M3S CIGAR strings will be grouped 
    /// if the coordinates of the matching substring (10M) 
    /// are the same between reads.
    #[arg(short='P', long="clip")]
    pub clip: bool,
    /// If enabled, reads will be grouped if their exon boundaries are the same. 
    /// This option discards any structural variants contained in mapped substrings of the read and only considers start and end coordinates of each non-splicing segment of the CIGAR string.
    #[arg(short='E', long="exon")]
    pub exon: bool,
    /// If enabled, supplementary alignments will be included in the collapsed groups of reads. 
    /// By default, TieBrush removes any mappings flagged as supplementary (0x800). Note, that if enabled, each supplementary mapping will count as a separate read
    #[arg(short='S', long="keep-supp")]
    pub keep_supplementary: bool,
    /// If enabled, secondary alignments will be included in the collapsed groups of reads. 
    /// By default, TieBrush removes any mappings not listed as primary (0x100). Note, that if enabled, each secondary mapping will count as a separate read
    #[arg(short='s', long="keep-secondary")]
    pub keep_secondary: bool,
    /// If enabled, unmapped reads will be included in the collapsed groups of reads. 
    /// By default, TieBrush removes any unmapped reads.
    #[arg(short='U', long="keep-unmapped")]
    pub keep_unmapped: bool,
    /// Maximum NH score of the reads to retain
    #[arg(short='N', long="max-nh")]
    pub max_nh: Option<u16>,
    /// Minimum quality score of the reads to retain
    #[arg(short='Q', long="min-qual")]
    pub min_qual: Option<u8>,
    /// Include only flags. 
    /// Reads with the provided bits set in the flag field are included.
    #[arg(short='f', long)]
    pub include_flags: Option<u16>,
}

trait ReadKeyStrat {
    type Key: Clone + PartialEq + Eq + std::hash::Hash + std::fmt::Debug;
    fn create_key(&self, record: &Record) -> anyhow::Result<Self::Key>;
}

#[derive(Debug, Copy, Clone)]
enum CmpStrat {
    Full,
    Clip,
    Exon,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct FullReadKey {
    tid: i32,
    start: i64,
    end: i64,
    strand: Strand,
    cigar: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ClipReadKey {
    tid: i32,
    start: i64,
    end: i64,
    strand: Strand,
    cigar: Vec<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ExonReadKey {
    tid: i32,
    strand: Strand,
    start: i64,
    end: i64,
    introns: Vec<(i64, i64)>,
}

struct FullStrat;
struct ClipStrat;
struct ExonStrat;

impl ReadKeyStrat for FullStrat {
    type Key = FullReadKey;

    fn create_key(&self, record: &Record) -> anyhow::Result<Self::Key> {
        let cigar = record.cigar().iter().map(|c| encode_cigar(c)).collect();
        Ok(Self::Key {
            tid: record.tid(),
            start: record.pos(),
            end: record.reference_end(),
            strand: get_strand(record)?,
            cigar,
        })
    }
}

impl ReadKeyStrat for ClipStrat {
    type Key = ClipReadKey;

    fn create_key(&self, record: &Record) -> anyhow::Result<Self::Key> {
        let cigar = record.cigar().iter()
            .filter(|c| !matches!(c, Cigar::SoftClip(_) | Cigar::HardClip(_)))
            .map(|c| encode_cigar(c))
            .collect();
        Ok(Self::Key {
            tid: record.tid(),
            start: record.pos(),
            end: record.reference_end(),
            strand: get_strand(record)?,
            cigar,
        })
    }
}

impl ReadKeyStrat for ExonStrat {
    type Key = ExonReadKey;

    fn create_key(&self, record: &Record) -> anyhow::Result<Self::Key> {
        let introns: Vec<(i64, i64)> = record.introns()
        .map(|arr|(arr[0], arr[1]))
        .collect();
        Ok(Self::Key {
            tid: record.tid(),
            strand: get_strand(record)?,
            start: record.pos(),
            end: record.reference_end(),
            introns,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum ReadKey {
    Full(FullReadKey),
    Clip(ClipReadKey),
    Exon(ExonReadKey),
}

impl ReadKey {
    fn create_key(tb_record: &TBSAMReaderRecord, cmp_strat: CmpStrat) -> anyhow::Result<Self> {
        match cmp_strat {
            CmpStrat::Full => FullStrat.create_key(tb_record.record()).map(ReadKey::Full),
            CmpStrat::Clip => ClipStrat.create_key(tb_record.record()).map(ReadKey::Clip),
            CmpStrat::Exon => ExonStrat.create_key(tb_record.record()).map(ReadKey::Exon),
        }
    }
}

#[derive(Debug)]
struct MergedReads {
    representative: Record,
    samples: Vec<bool>,
    acc_yc: u32,
    acc_yx: u32,
}

impl MergedReads {
    fn new(tb_record: &TBSAMReaderRecord, num_samples: usize) -> Self {
        let yc = get_yc_tag(tb_record.record()).expect("Failed to read YC tag").unwrap_or(1) as u32;
        let yx = get_yx_tag(tb_record.record()).expect("Failed to read YX tag").unwrap_or(1) as u32;
        let sample_idx = tb_record.file_idx();
        let mut samples = vec![false; num_samples];
        samples[sample_idx] = true;
        Self {
            representative: tb_record.record().clone(),
            samples,
            acc_yc: yc,
            acc_yx: yx,
        }
    }

    fn add_duplicate(&mut self, tb_record: &TBSAMReaderRecord) {
        let yc = get_yc_tag(tb_record.record()).expect("Failed to read YC tag").unwrap_or(1) as u32;
        let yx = get_yx_tag(tb_record.record()).expect("Failed to read YX tag").unwrap_or(1) as u32;
        let sample_idx = tb_record.file_idx();

        self.samples[sample_idx] = true;
        self.acc_yc += yc;
        self.acc_yx += yx;
    }
}

pub struct BrushCMD {
    brush_args: BrushArgs,
    tb_writer: Option<Writer>,
    tb_format: Format,
    cmp_strat: CmpStrat,

    current_reads: HashMap<ReadKey, MergedReads>,
    last_position: Option<(i32, i64)>, // (tid, pos)

    tb_stats: TBStats,
}

impl BrushCMD {
    pub fn new(args: BrushArgs) -> anyhow::Result<Self> {
        // verify all provided alignments exists
        for alignment in &args.input_alignments {
            if !alignment.exists() {
                anyhow::bail!("Alignment file {} does not exist", alignment.display());
            }
        }

        // get output format
        let tb_format = get_format(args.output.as_ref().unwrap())?;
        match tb_format {
            Format::Sam | Format::Bam => {}
            _ => anyhow::bail!("Unsupported output format. Make sure the output file has a .sam or .bam extension"),
        }

        let cmp_strat = if args.full {
            CmpStrat::Full
        } else if args.clip {
            CmpStrat::Clip
        } else if args.exon {
            CmpStrat::Exon
        } else {
            CmpStrat::Full
        };

        Ok(Self {
            brush_args: args,
            tb_writer: None,
            tb_format,
            cmp_strat,
            current_reads: HashMap::new(),
            last_position: None,
            tb_stats: TBStats::default(),
        })
    }

    pub fn run(&mut self) -> anyhow::Result<()> {
        let req_fields = htslib::sam_fields_SAM_QNAME | htslib::sam_fields_SAM_FLAG | htslib::sam_fields_SAM_RNAME | htslib::sam_fields_SAM_POS | htslib::sam_fields_SAM_CIGAR | htslib::sam_fields_SAM_AUX;
        let sam_reader = TBSAMReader::new_with_fields(&self.brush_args.input_alignments, req_fields)?;
        let header = sam_reader.get_header();

        // write TB header out
        let tb_path = self.brush_args.output.as_ref().unwrap();
        let writer = Writer::from_path(tb_path, header, self.tb_format)?;
        self.tb_writer = Some(writer);

        for tb_record in sam_reader {
            self.tb_stats.total_input_reads += 1;
            // check if the record passes the options
            if !self.brush_args.keep_supplementary && tb_record.record().is_supplementary() {
                self.tb_stats.filtered_reads += 1;
                continue;
            }
            if !self.brush_args.keep_secondary && tb_record.record().is_secondary() {
                self.tb_stats.filtered_reads += 1;
                continue;
            }
            if !self.brush_args.keep_unmapped && tb_record.record().is_unmapped() {
                self.tb_stats.filtered_reads += 1;
                continue;
            }
            if tb_record.record().mapq() < self.brush_args.min_qual.unwrap_or(0) {
                self.tb_stats.filtered_reads += 1;
                continue;
            }
            if get_nh_tag(&tb_record.record())?.unwrap_or(1) > self.brush_args.max_nh.unwrap_or(u16::MAX) {
                self.tb_stats.filtered_reads += 1;
                continue;
            }
            if self.brush_args.include_flags.is_some() && !flags_set(&tb_record.record(), self.brush_args.include_flags.unwrap()) {
                self.tb_writer.as_mut().unwrap().write(&tb_record.record())?;
                self.tb_stats.total_output_reads += 1;
                continue;
            }

            // everything else can be processed via merging
            self.tb_stats.processed_reads += 1;
            self.process_record(tb_record)?;

        }
        self.flush_current_reads()?;

        // Print summary
        self.tb_stats.print_summary();

        Ok(())
    }

    fn process_record(&mut self, record: TBSAMReaderRecord) -> anyhow::Result<()> {
        let current_pos = (record.record().tid(), record.record().pos());
        
        match self.last_position {
            Some((tid, pos)) if tid == current_pos.0 && pos == current_pos.1 => {},
            _ => {
                self.flush_current_reads()?;
            }
        }
        
        self.last_position = Some(current_pos);
        
        let read_key = ReadKey::create_key(&record,self.cmp_strat)?;
        if let Some(grouped_reads) = self.current_reads.get_mut(&read_key) {
            grouped_reads.add_duplicate(&record);
        } else {
            self.current_reads.insert(read_key, MergedReads::new(&record, self.brush_args.input_alignments.len()));
        }
        
        Ok(())
    }

    fn flush_current_reads(&mut self) -> anyhow::Result<()> {
        for (_, mut grouped_reads) in self.current_reads.drain() {
            grouped_reads.representative.push_aux(b"YC", rust_htslib::bam::record::Aux::I32(grouped_reads.acc_yc as i32))?;
            grouped_reads.representative.push_aux(b"YX",rust_htslib::bam::record::Aux::I32(grouped_reads.samples.iter().filter(|&&s| s).count() as i32))?;
            
            self.tb_writer.as_mut().unwrap().write(&grouped_reads.representative)?;
            self.tb_stats.total_output_reads += 1;
        }
        Ok(())
    }
}

pub fn run(args: BrushArgs) -> anyhow::Result<()> {
    let mut cmd = BrushCMD::new(args)?;
    cmd.run()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{Record, Header, header::HeaderRecord};
    use rust_htslib::bam::record::{Cigar, CigarString};

    fn create_test_header() -> Header {
        let mut header = Header::new();
        
        // Add HD record first (required)
        header.push_record(&HeaderRecord::new(b"HD")
            .push_tag(b"VN", "1.6")
            .push_tag(b"SO", "coordinate"));
        
        // Add SQ record for reference sequence
        header.push_record(&HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", "1000000"));  // LN should be bytes, not &str
        
        header
    }

    fn create_test_record(read_name: &[u8], tid: i32, pos: i64, cigar: Vec<Cigar>, strand: bool) -> Record {
        let mut record = Record::new();
        
        record.set_tid(tid);
        record.set_pos(pos);
        
        let seq_len = 100;
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        record.set(read_name, Some(&CigarString(cigar)), &seq, &qual);
        
        record.set_mapq(60);
        
        let mut flags = 0u16;
        if !strand {
            flags |= 16;
        }
        record.set_flags(flags);
        
        record.set_mtid(-1);
        record.set_mpos(-1);
        record.set_insert_size(0);
        
        record
    }

    #[test]
    fn test_full_strategy_same_reads() {
        let strat = FullStrat;
        let cigar = vec![
            Cigar::Match(50), 
            Cigar::SoftClip(10)
        ];
        
        let record1 = create_test_record(b"r1",0, 100, cigar.clone(), true);
        let record2 = create_test_record(b"r2",0, 100, cigar.clone(), true);
        
        let key1 = strat.create_key(&record1).unwrap();
        let key2 = strat.create_key(&record2).unwrap();
        
        assert_eq!(key1, key2);
    }

    #[test]
    fn test_full_strategy_different_positions() {
        let strat = FullStrat;
        let cigar = vec![
            Cigar::Match(50), 
            Cigar::SoftClip(10)
        ];
        
        let record1 = create_test_record(b"r1",0, 100, cigar.clone(), true);
        let record2 = create_test_record(b"r2",0, 200, cigar.clone(), true);
        
        let key1 = strat.create_key(&record1).unwrap();
        let key2 = strat.create_key(&record2).unwrap();
        
        assert_ne!(key1, key2);
    }

    #[test]
    fn test_merged_reads_fresh_records() {
        // Create mock TBSAMReaderRecord
        let record1 = create_test_record(b"r1",0, 100, vec![Cigar::Match(50)], true);
        let record2 = create_test_record(b"r2",0, 100, vec![Cigar::Match(50)], true);
    
        let tb_rec1 = TBSAMReaderRecord::new(record1, 0);
        let tb_rec2 = TBSAMReaderRecord::new(record2, 1);

        let num_samples = 2;

        let mut merged_reads = MergedReads::new(&tb_rec1, num_samples);
        assert_eq!(merged_reads.samples.len(), num_samples);
        assert!(merged_reads.samples[0]);
        assert!(!merged_reads.samples[1]);
        assert_eq!(merged_reads.acc_yc, 1);
        assert_eq!(merged_reads.acc_yx, 1);
        
        merged_reads.add_duplicate(&tb_rec2);
        assert_eq!(merged_reads.samples.len(), num_samples);
        assert!(merged_reads.samples[0]);
        assert!(merged_reads.samples[1]);
        assert_eq!(merged_reads.acc_yc, 2);
        assert_eq!(merged_reads.acc_yx, 2);
    }

    #[test]
    fn test_merged_reads_with_existing_tags() {
        let record1 = create_test_record(b"r1", 0, 100, vec![Cigar::Match(50)], true);
        let mut record2 = create_test_record(b"r2", 0, 100, vec![Cigar::Match(50)], true);
        record2.push_aux(b"YC", rust_htslib::bam::record::Aux::I32(4)).unwrap();
        record2.push_aux(b"YX", rust_htslib::bam::record::Aux::I32(2)).unwrap();

        let tb_rec1 = TBSAMReaderRecord::new(record1, 0);
        let tb_rec2 = TBSAMReaderRecord::new(record2, 1);

        let num_samples = 2;

        let mut merged_reads = MergedReads::new(&tb_rec1, num_samples);
        assert_eq!(merged_reads.samples.len(), num_samples);
        assert!(merged_reads.samples[0]);
        assert!(!merged_reads.samples[1]);
        assert_eq!(merged_reads.acc_yc, 1);
        assert_eq!(merged_reads.acc_yx, 1);

        merged_reads.add_duplicate(&tb_rec2);
        assert_eq!(merged_reads.samples.len(), num_samples);
        assert!(merged_reads.samples[0]);
        assert!(merged_reads.samples[1]);
        assert_eq!(merged_reads.acc_yc, 5);
        assert_eq!(merged_reads.acc_yx, 3);
    }
}