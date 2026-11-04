//! Recording this tool's invocation as an `@PG` header line, as samtools does.

use std::env;

use anyhow::Context;
use pasta::sam::{
    self as sam,
    header::record::value::{
        Map,
        map::{Program, program::tag},
    },
};

const NAME: &str = "ntools";

/// Adds an `@PG` line for this invocation at the end of each program chain in the header, as
/// samtools does (`sam_hdr_add_pg`).
pub fn add(header: &mut sam::Header) -> anyhow::Result<()> {
    let command_line = env::args().collect::<Vec<_>>().join(" ");

    let map = Map::<Program>::builder()
        .insert(tag::NAME, NAME)
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(tag::COMMAND_LINE, command_line)
        .build()
        .context("building @PG record")?;

    header
        .add_program_record(NAME, map)
        .context("adding @PG record")
}
