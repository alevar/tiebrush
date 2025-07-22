use crate::samreader::SAMReader;

use clap::{Args, ArgGroup};
use std::path::PathBuf;

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

pub fn run(args: CovArgs) -> anyhow::Result<()> {
    // verify all provided alignments exists
    for alignment in &args.input_alignments {
        if !alignment.exists() {
            anyhow::bail!("Alignment file {} does not exist", alignment.display());
        }
    }

    // verify that coverage exists
    match args.coverage {
        Some(coverage) => {
            if !coverage.exists() {
                anyhow::bail!("Coverage file {} does not exist", coverage.display());
            }
        }
        None => {}
    }

    // verify that junctions exists
    match args.junctions {
        Some(junctions) => {
            if !junctions.exists() {
                anyhow::bail!("Junctions file {} does not exist", junctions.display());
            }
        }
        None => {}
    }

    // read inputs
    let sam_reader = SAMReader::new(&args.input_alignments)?;
    // debug by printing the first N records
    for (i, record) in sam_reader.take(1000).enumerate() {
        println!("{:?}", record);
    }
    
    Ok(())
}