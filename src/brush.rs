use crate::samreader::TBSAMReader;
use crate::commons::*;

use std::{path::PathBuf, vec};
use std::collections::HashMap;
use clap::{Args, ArgGroup};
use rust_htslib::bam::{Record, record::Cigar, Format, ext::BamRecordExtensions, Writer};

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
    /// By default, TieBrush removes any mappings not listed as primary (0x100). Note, that if enabled, each supplementary mapping will count as a separate read
    #[arg(short='S', long="keep-supp")]
    pub keep_supplementary: bool,
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
    segments: Vec<(i64, i64)>,
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
            end: record.cigar().end_pos(),
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
            end: record.cigar().end_pos(),
            strand: get_strand(record)?,
            cigar,
        })
    }
}

impl ReadKeyStrat for ExonStrat {
    type Key = ExonReadKey;

    fn create_key(&self, record: &Record) -> anyhow::Result<Self::Key> {
        let segments = vec![];
        Ok(Self::Key {
            tid: record.tid(),
            strand: get_strand(record)?,
            segments,
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
    fn create_key(record: &Record, cmp_strat: CmpStrat) -> anyhow::Result<Self> {
        match cmp_strat {
            CmpStrat::Full => FullStrat.create_key(record).map(ReadKey::Full),
            CmpStrat::Clip => ClipStrat.create_key(record).map(ReadKey::Clip),
            CmpStrat::Exon => ExonStrat.create_key(record).map(ReadKey::Exon),
        }
    }
}

#[derive(Debug)]
struct MergedReads {
    representative: Record,
    dup_count: u32,
    sample_count: u32,
}

impl MergedReads {
    fn new(record: Record) -> Self {
        Self {
            representative: record,
            dup_count: 1,
            sample_count: 1,
        }
    }

    fn add_duplicate(&mut self) {
        self.dup_count += 1;
    }
}

pub struct BrushCMD {
    brush_args: BrushArgs,
    tb_writer: Option<Writer>,
    tb_format: Format,
    cmp_strat: CmpStrat,

    current_reads: HashMap<ReadKey, MergedReads>,
    last_position: Option<(i32, i64)>, // (tid, pos)
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
        })
    }

    pub fn run(&mut self) -> anyhow::Result<()> {
        let sam_reader = TBSAMReader::new(&self.brush_args.input_alignments)?;
        let header = sam_reader.get_header();

        // write TB header out
        let tb_path = self.brush_args.output.as_ref().unwrap();
        let writer = Writer::from_path(tb_path, header, self.tb_format)?;
        self.tb_writer = Some(writer);

        for record in sam_reader {
            // check if the record passes the options
            if !self.brush_args.keep_supplementary && record.is_supplementary() {
                continue;
            }
            if !self.brush_args.keep_unmapped && record.is_unmapped() {
                continue;
            }
            if record.mapq() < self.brush_args.min_qual.unwrap_or(0) {
                continue;
            }
            if get_nh_tag(&record)?.unwrap_or(1) > self.brush_args.max_nh.unwrap_or(u16::MAX) {
                continue;
            }
            if self.brush_args.include_flags.is_some() && !flags_set(&record, self.brush_args.include_flags.unwrap()) {
                self.tb_writer.as_mut().unwrap().write(&record)?;
                continue;
            }

            // everything else can be processed via merging
            self.process_record(record)?;

        }
        self.flush_current_reads()?;

        Ok(())
    }

    fn process_record(&mut self, record: Record) -> anyhow::Result<()> {
        let current_pos = (record.tid(), record.pos());
        
        match self.last_position {
            Some((tid, pos)) if tid == current_pos.0 && pos == current_pos.1 => {},
            _ => {
                self.flush_current_reads()?;
            }
        }
        
        self.last_position = Some(current_pos);
        
        let read_key = ReadKey::create_key(&record,self.cmp_strat)?;
        if let Some(grouped_reads) = self.current_reads.get_mut(&read_key) {
            grouped_reads.add_duplicate();
        } else {
            self.current_reads.insert(read_key, MergedReads::new(record));
        }
        
        Ok(())
    }

    fn flush_current_reads(&mut self) -> anyhow::Result<()> {
        for (_, mut grouped_reads) in self.current_reads.drain() {
            if grouped_reads.dup_count > 1 {
                grouped_reads.representative.push_aux(b"YC", rust_htslib::bam::record::Aux::I32(grouped_reads.dup_count as i32))?;
            }
            
            self.tb_writer.as_mut().unwrap().write(&grouped_reads.representative)?;
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

    fn create_test_record(tid: i32, pos: i64, cigar: Vec<Cigar>, strand: bool) -> Record {
        let header = create_test_header();
        let mut record = Record::new();
        
        record.set_tid(tid);
        record.set_pos(pos);
        
        let seq_len = 100;
        let seq = vec![b'A'; seq_len];
        let qual = vec![30u8; seq_len];
        record.set(b"test_read", Some(&CigarString(cigar)), &seq, &qual);
        
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
        let cigar = vec![Cigar::Match(50), Cigar::SoftClip(10)];
        
        let record1 = create_test_record(0, 100, cigar.clone(), true);
        let record2 = create_test_record(0, 100, cigar.clone(), true);
        
        let key1 = strat.create_key(&record1).unwrap();
        let key2 = strat.create_key(&record2).unwrap();
        
        assert_eq!(key1, key2);
    }

    #[test]
    fn test_full_strategy_different_positions() {
        let strat = FullStrat;
        let cigar = vec![Cigar::Match(50), Cigar::SoftClip(10)];
        
        let record1 = create_test_record(0, 100, cigar.clone(), true);
        let record2 = create_test_record(0, 200, cigar.clone(), true);
        
        let key1 = strat.create_key(&record1).unwrap();
        let key2 = strat.create_key(&record2).unwrap();
        
        assert_ne!(key1, key2);
    }
}