use crate::samreader::TBSAMReader;
use crate::commons::*;

use std::fs::File;
use std::path::PathBuf;
use clap::{Args, ArgGroup};
use std::collections::{HashMap, hash_map::Entry};
use rust_htslib::bam::{Record, HeaderView, record::Cigar, ext::BamRecordExtensions};

#[derive(Args, Debug)]
pub struct TestArgs {
    /// One or more alignment files in SAM/BAM/CRAM format.
    #[arg(required = true)]
    pub input_alignments: Vec<PathBuf>,
}

pub struct TestCMD {
    test_args: TestArgs,
    sam_reader: TBSAMReader,
    header_view: HeaderView,
}

impl TestCMD {
    pub fn new(args: TestArgs) -> anyhow::Result<Self> {
        // verify all provided alignments exists
        for alignment in &args.input_alignments {
            if !alignment.exists() {
                anyhow::bail!("Alignment file {} does not exist", alignment.display());
            }
        }

        // Create SAM reader to get header for BigWig initialization
        let sam_reader = TBSAMReader::new(&args.input_alignments)?;
        let header = sam_reader.get_header();
        let header_view = HeaderView::from_header(header);

        Ok(Self {
            test_args: args,
            sam_reader,
            header_view,
        })
    }

    pub fn run(&mut self) -> anyhow::Result<()> {
        let mut record_count = 0;
        while let Some(tb_record) = self.sam_reader.next() {
            let record = tb_record.record();
            record_count += 1;
        }
        println!("Total records processed: {}", record_count);

        Ok(())
    }
}

pub fn run(args: TestArgs) -> anyhow::Result<()> {
    let mut cmd = TestCMD::new(args)?;
    cmd.run()?;

    Ok(())
}