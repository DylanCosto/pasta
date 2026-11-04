//! Statistics of two inputs (`bcftools stats A B`): the synced reader pairs their records,
//! which are counted in three sets (only in A, only in B, in both), and with `-s/-S` the
//! genotypes of the samples both inputs have are compared at the paired records.

use std::{
    io::{self, Write},
    process::ExitCode,
};

use anyhow::bail;
use pasta::variant::{
    Unpack,
    synced::{IndexRequirement, Pairing, SyncedOptions, SyncedReaders, TargetsLogic},
};

use super::{
    Args, OUT_BUFFER, Run, Samples, Setup, check_options,
    collect::{Collector, Fatal},
    compare::{GtCompare, Side},
    compile_filter, fatal_exit,
    indel_context::IndelContext,
    init_stats, open_error, open_exons, read_error,
    report::{self, ReportInfo},
};
use crate::filter::FilterLogic;

/// The pairing of `-c/--collapse` values, in order (bcftools' `files->collapse`).
fn pairing(collapse: &[String]) -> Pairing {
    let mut pair = Pairing::empty();
    for c in collapse {
        match c.as_str() {
            "snps" => pair |= Pairing::SNPS,
            "indels" => pair |= Pairing::INDELS,
            "both" => pair |= Pairing::BOTH,
            "any" | "all" => pair |= Pairing::ANY,
            "some" => pair |= Pairing::SOME,
            _ => pair = Pairing::empty(),
        }
    }
    pair
}

/// Opens both inputs with the options of the command line (bcftools' `main_vcfstats`): they
/// must be indexed.
fn open(args: &Args, files: &[String], options: SyncedOptions) -> anyhow::Result<SyncedReaders> {
    let mut sr = SyncedReaders::new();
    sr.set_index_requirement(IndexRequirement::Required);
    sr.set_pairing(pairing(&args.collapse));
    if let Some(t) = options.targets {
        sr.set_targets_overlap(t.overlap);
        let logic = if options.targets_exclude {
            TargetsLogic::Exclude
        } else {
            TargetsLogic::Include
        };
        sr.set_targets(t, logic, None)?;
    }
    if let Some(r) = options.regions {
        sr.set_regions_overlap(r.overlap);
        sr.set_regions(r);
    }
    sr.set_threads(options.threads);
    sr.set_apply_filters(options.apply_filters.as_deref());
    for f in files {
        if let Err(e) = sr.add_reader(f) {
            return Err(open_error(f, &e));
        }
    }
    Ok(sr)
}

/// Resolves `-s/-S` to the samples both inputs have (`bcf_sr_set_samples`).
fn select_samples(sr: &mut SyncedReaders, run: &Run<'_>) -> anyhow::Result<Samples> {
    let Some((spec, source)) = run.sample_spec else {
        return Ok(Samples {
            columns: vec![Vec::new(); sr.reader_count()],
            names: Vec::new(),
            given: false,
        });
    };
    // A list that cannot be read was reported by the synced reader; it fails as an empty one.
    if !sr.set_samples(spec.as_bytes(), source).unwrap_or(false) {
        if sr.header(0).sample_count() == 0 {
            bail!("No sample columns in {}", sr.file_name(0));
        }
        bail!("Unable to parse the samples: \"{spec}\"");
    }
    let columns = (0..sr.reader_count())
        .map(|i| sr.sample_map(i).iter().map(|&c| c as usize).collect())
        .collect();
    Ok(Samples {
        columns,
        names: sr.samples().to_vec(),
        given: true,
    })
}

/// Runs `stats A B`.
pub(super) fn run(args: &Args, run: Run<'_>, options: SyncedOptions) -> anyhow::Result<ExitCode> {
    let files = run.files;
    let verbose_sites = run.verbose_sites;
    let mut sr = open(args, files, options)?;
    let bins = check_options(args, &[sr.header(0), sr.header(1)], run.filter.as_ref())?;
    let samples = select_samples(&mut sr, &run)?;
    let setup = init_stats(args, &[sr.header(0), sr.header(1)], samples, bins, run)?;
    let exons = open_exons(args)?;
    sr.set_max_unpack(setup.max_unpack());
    let indel_ctx = setup.fasta_ref.as_deref().and_then(IndelContext::open);

    let nsamples = [sr.header(0).sample_count(), sr.header(1).sample_count()];
    let mut out = io::BufWriter::with_capacity(OUT_BUFFER, crate::output::stdout());
    let info = ReportInfo {
        files,
        nsamples: &nsamples,
        sample_names: &setup.sample_names,
        split_by_id: false,
        verbose_sites,
        exons: args.exons.is_some(),
        indel_ctx: indel_ctx.is_some(),
    };
    report::print_header(&mut out, &info)?;

    let cfg = &setup.cfg;
    let mut collector = match Collector::try_new(cfg, &setup.user_types, indel_ctx, exons) {
        Ok(collector) => collector,
        Err(f) => return fatal_exit(&mut out, f),
    };
    let mut gtcmp = (cfg.n_smpl() > 0).then(|| GtCompare::new(cfg.m_af, cfg.n_smpl()));
    let verbose = verbose_sites.then_some(&mut out);
    let collected = collect(&mut sr, &setup, &mut collector, gtcmp.as_mut(), verbose);
    if let Err(f) = collected {
        return fatal_exit(&mut out, f);
    }
    if let Some(e) = sr.take_pending_error() {
        out.flush()?;
        bail!("Error: {e}");
    }
    report::print_stats(&mut out, cfg, &info, &mut collector.sets, gtcmp.as_mut())?;
    out.flush()?;
    Ok(ExitCode::SUCCESS)
}

/// Collects the statistics of every position (`do_vcf_stats` with two readers): a record
/// counts for the set of the inputs that have it, unless the filter rejects it in any input.
/// With `verbose`, the genotype comparison of each site is written to it.
fn collect(
    sr: &mut SyncedReaders,
    setup: &Setup,
    collector: &mut Collector,
    mut gtcmp: Option<&mut GtCompare>,
    mut verbose: Option<&mut impl Write>,
) -> Result<(), Fatal> {
    let cfg = &setup.cfg;
    let mut filters = [
        compile_filter(setup, sr.header(0))?,
        compile_filter(setup, sr.header(1))?,
    ];
    let write_error = |e: io::Error| Fatal::error(format!("Error: {e}\n"));
    loop {
        match sr.next_line() {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(read_error(&e)),
        }
        let mut matched = 0usize;
        let mut first = None;
        let mut pass = true;
        for (i, filter) in filters.iter_mut().enumerate() {
            let Some((rec, header)) = sr.line_with_header(i) else {
                continue;
            };
            if let Some((f, logic)) = filter {
                let ok = f
                    .pass(header, rec)
                    .map_err(|e| Fatal::error(format!("{e:#}\n")))?;
                if ok != (*logic == FilterLogic::Include) {
                    pass = false;
                    break;
                }
            }
            matched |= 1 << i;
            first.get_or_insert(i);
        }
        let Some(reader) = first.filter(|_| pass) else {
            continue;
        };
        let (rec, header) = sr.line_with_header(reader).expect("the reader has a line");
        collector.add(cfg, header, rec, matched - 1, reader)?;

        // Genotypes are compared at the first ALT allele's frequency, which a record without
        // ALT alleles keeps from an earlier record, as bcftools does.
        if matched == 3
            && let Some(g) = gtcmp.as_deref_mut()
            && let Some(iaf) = collector.first_alt_bin()
        {
            let line_type = sr.line_mut(0).expect("paired").variant_types();
            sr.line_mut(1).expect("paired").unpack(Unpack::FMT);
            let a = Side {
                header: sr.header(0),
                rec: sr.line(0).expect("paired"),
                samples: &cfg.samples[0],
            };
            let b = Side {
                header: sr.header(1),
                rec: sr.line(1).expect("paired"),
                samples: &cfg.samples[1],
            };
            let verbose = verbose
                .as_deref_mut()
                .map(|out| (out as &mut dyn Write, &setup.sample_names[..]));
            g.add(&a, &b, iaf, line_type, verbose)
                .map_err(write_error)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collapse_values_accumulate_until_none() {
        let p = |v: &[&str]| pairing(&v.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>());
        assert_eq!(p(&[]), Pairing::empty());
        assert_eq!(p(&["snps", "indels"]), Pairing::BOTH);
        assert_eq!(p(&["all"]), Pairing::ANY);
        assert_eq!(p(&["some", "none", "snps"]), Pairing::SNPS);
    }
}
