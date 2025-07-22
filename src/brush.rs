use clap::Parser;
use std::path::PathBuf;

#[derive(Parser)]
pub struct BrushArgs {
    #[arg(short, long)]
    pub input: PathBuf,
}