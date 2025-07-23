use crate::samreader::TBSAMReader;
use crate::commons::*;

use std::fs::File;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use clap::{Args, ArgGroup};
use std::collections::{HashMap, hash_map::Entry};
use rust_htslib::bam::{Record, HeaderView, record::Cigar, Format, ext::BamRecordExtensions, Writer};

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

pub struct BrushCMD {
    brush_args: BrushArgs,
    tb_writer: Option<Writer>,
    tb_format: Format,
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

        Ok(Self {
            brush_args: args,
            tb_writer: None,
            tb_format,
        })
    }

    pub fn run(&mut self) -> anyhow::Result<()> {
        let mut sam_reader = TBSAMReader::new(&self.brush_args.input_alignments)?;
        let header = sam_reader.get_header();

        // write TB header out
        let tb_path = self.brush_args.output.as_ref().unwrap();
        let writer = Writer::from_path(tb_path, header, self.tb_format)?;
        self.tb_writer = Some(writer);

        Ok(())
    }
}

pub fn run(args: BrushArgs) -> anyhow::Result<()> {
    let mut cmd = BrushCMD::new(args)?;
    cmd.run()?;

    Ok(())
}