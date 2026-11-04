//! `+counts`: the number of samples, SNPs, indels, MNPs, other variants, and sites (bcftools'
//! `plugins/counts.c`).

use std::io::Write;

use pasta::variant::{Header, Record, Unpack, VariantType};

use super::{InitHeaders, Loaded, Output, Plugin};

/// The description `plugin -lv` prints.
pub const ABOUT: &str = "A minimal plugin which counts number of samples, SNPs,\nINDELs, MNPs and total number of sites.\n";

/// The usage text (bcftools prints the description, which is all the plugin has).
pub fn usage() -> String {
    ABOUT.to_owned()
}

/// The counts (bcftools' globals).
#[derive(Debug, Default)]
struct Counts {
    nsamples: u64,
    nsnps: u64,
    nindels: u64,
    nmnps: u64,
    nothers: u64,
    nsites: u64,
}

/// Starts counting; nothing is written but the counts (bcftools' `init` returning 1).
pub fn init(_args: &[String], hdrs: InitHeaders<'_>) -> anyhow::Result<Loaded> {
    let Some((hdr, _)) = hdrs else {
        // bcftools reads the sample count of a missing header and crashes.
        return Err(super::general_usage());
    };
    Ok(Loaded {
        plugin: Box::new(Counts {
            nsamples: hdr.sample_count() as u64,
            ..Counts::default()
        }),
        output: Output::None,
        // Only the alleles are read: the VCF parser skips the rest of each line, where
        // bcftools parses it all.
        max_unpack: Unpack::STR,
    })
}

impl Plugin for Counts {
    fn process(&mut self, _hdr: &Header, _out: &Header, rec: &mut Record) -> anyhow::Result<bool> {
        rec.unpack(Unpack::STR);
        let ty = rec.variant_types();
        self.nsnps += u64::from(ty.contains(VariantType::SNP));
        self.nindels += u64::from(ty.contains(VariantType::INDEL));
        self.nmnps += u64::from(ty.contains(VariantType::MNP));
        self.nothers += u64::from(ty.contains(VariantType::OTHER));
        self.nsites += 1;
        Ok(false)
    }

    fn destroy(&mut self) -> anyhow::Result<()> {
        let mut out = crate::output::stdout();
        write!(
            out,
            "Number of samples: {}\nNumber of SNPs:    {}\nNumber of INDELs:  {}\nNumber of MNPs:    {}\nNumber of others:  {}\nNumber of sites:   {}\n",
            self.nsamples, self.nsnps, self.nindels, self.nmnps, self.nothers, self.nsites
        )?;
        out.flush()?;
        Ok(())
    }
}
