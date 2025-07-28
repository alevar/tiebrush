mod brush;
mod cov;
#[cfg(feature = "test")]
mod test;
mod samreader;
mod commons;

use clap::{Parser, Subcommand};
use anyhow;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(subcommand_required = true, arg_required_else_help = true)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Brush(brush::BrushArgs),
    Cov(cov::CovArgs),
    #[cfg(feature = "test")]
    Test(test::TestArgs),
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Brush(brush_ops ) => {
            brush::run(brush_ops)?;
        }
        Commands::Cov(cov_ops) => {
            cov::run(cov_ops)?;
        }
        #[cfg(feature = "test")]
        Commands::Test(test_ops) => {
            test::run(test_ops)?;
        }
    }

    Ok(())
}