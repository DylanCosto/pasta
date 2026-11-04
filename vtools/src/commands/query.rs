//! `vtools query`: print fields of VCF/BCF records in a user-defined format (bcftools query).

mod inputs;

use std::{
    fs::File,
    io::{self, BufWriter, IsTerminal, Write},
    path::PathBuf,
    process::ExitCode,
};

use anyhow::{Context, bail};
use clap::ArgMatches;
use pasta::variant::{
    Header, Record, Unpack,
    readlist::{self, ListSource},
    synced::{Overlap, Regions, SyncedOptions, SyncedReader},
};

use self::inputs::Inputs;

use crate::options::{last_index, later};
use crate::{
    commands::head::apply_verbosity,
    convert::Convert,
    filter::{self, FilterLogic, RecordFilter},
    raw_input::strerror,
};

#[derive(Debug, clap::Args)]
#[command(disable_help_flag = true, args_override_self = true)]
pub struct Args {
    /// Input VCF/BCF files (`-` or none for stdin); with several, records of the first are
    /// printed and all need an index.
    inputs: Vec<String>,

    /// Exclude sites for which the expression is true.
    #[arg(
        short = 'e',
        long = "exclude",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    exclude: Vec<String>,

    /// Select sites for which the expression is true.
    #[arg(
        short = 'i',
        long = "include",
        value_name = "EXPR",
        allow_hyphen_values = true
    )]
    include: Vec<String>,

    /// Only warn about unknown subset samples.
    #[arg(long = "force-samples")]
    force_samples: bool,

    /// Output STR for samples failing the -i/-e filtering expression.
    #[arg(
        short = 'F',
        long = "print-filtered",
        value_name = "STR",
        allow_hyphen_values = true
    )]
    print_filtered: Option<String>,

    /// The output format, see bcftools' documentation.
    #[arg(
        short = 'f',
        long = "format",
        value_name = "STRING",
        allow_hyphen_values = true
    )]
    format: Option<String>,

    /// Print header, -HH to omit column indices.
    #[arg(short = 'H', long = "print-header", action = clap::ArgAction::Count)]
    print_header: u8,

    /// Print the list of samples and exit.
    #[arg(short = 'l', long = "list-samples")]
    list_samples: bool,

    /// Disable automatic addition of newline character when not present.
    #[arg(short = 'N', long = "disable-automatic-newline")]
    no_newline: bool,

    /// Output file name [stdout].
    #[arg(
        short = 'o',
        long = "output",
        visible_alias = "output-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    output: Option<PathBuf>,

    /// Restrict to comma-separated list of regions.
    #[arg(
        short = 'r',
        long = "regions",
        value_name = "REGION",
        allow_hyphen_values = true
    )]
    regions: Option<String>,

    /// Restrict to regions listed in a file.
    #[arg(
        short = 'R',
        long = "regions-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    regions_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1].
    #[arg(
        long = "regions-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true
    )]
    regions_overlap: Option<String>,

    /// Similar to -r but streams rather than index-jumps.
    #[arg(
        short = 't',
        long = "targets",
        value_name = "[^]REGION",
        allow_hyphen_values = true
    )]
    targets: Option<String>,

    /// Similar to -R but streams rather than index-jumps.
    #[arg(
        short = 'T',
        long = "targets-file",
        value_name = "[^]FILE",
        allow_hyphen_values = true
    )]
    targets_file: Option<String>,

    /// Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0].
    #[arg(
        long = "targets-overlap",
        value_name = "0|1|2",
        allow_hyphen_values = true
    )]
    targets_overlap: Option<String>,

    /// Print CHROM, POS, MASK, REF, ALT and the listed INFO tags (comma-separated).
    #[arg(
        short = 'a',
        long = "annots",
        value_name = "LIST",
        hide = true,
        allow_hyphen_values = true
    )]
    annots: Option<String>,

    /// List of samples to include.
    #[arg(
        short = 's',
        long = "samples",
        value_name = "LIST",
        allow_hyphen_values = true
    )]
    samples: Option<String>,

    /// File of samples to include.
    #[arg(
        short = 'S',
        long = "samples-file",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    samples_file: Option<String>,

    /// Obsolete.
    #[arg(
        short = 'c',
        long = "collapse",
        value_name = "STRING",
        hide = true,
        allow_hyphen_values = true
    )]
    collapse: Option<String>,

    /// Process multiple VCFs listed in the file.
    #[arg(
        short = 'v',
        long = "vcf-list",
        value_name = "FILE",
        allow_hyphen_values = true
    )]
    vcf_list: Option<String>,

    /// Print "." for undefined tags.
    #[arg(short = 'u', long = "allow-undef-tags")]
    allow_undef_tags: bool,

    /// Verbosity level.
    #[arg(long = "verbosity", value_name = "INT", allow_hyphen_values = true)]
    verbosity: Option<String>,

    /// Use multithreading with INT worker threads for decompression [0] (a vtools extension).
    #[arg(long = "threads", value_name = "INT", allow_hyphen_values = true)]
    threads: Option<String>,

    /// Print help.
    #[arg(short = 'h', long = "help")]
    help: bool,
}

const USAGE: &str = "
About:   Extracts fields from VCF/BCF file and prints them in user-defined format
Usage:   vtools query [options] <A.vcf.gz> [<B.vcf.gz> [...]]

Options:
    -e, --exclude EXPR                Exclude sites for which the expression is true (see man page for details)
        --force-samples               Only warn about unknown subset samples
    -F, --print-filtered STR          Output STR for samples failing the -i/-e filtering expression
    -f, --format STRING               See man page for details
    -H, --print-header                Print header, -HH to omit column indices
    -i, --include EXPR                Select sites for which the expression is true (see man page for details)
    -l, --list-samples                Print the list of samples and exit
    -N, --disable-automatic-newline   Disable automatic addition of newline character when not present
    -o, --output FILE                 Output file name [stdout]
    -r, --regions REGION              Restrict to comma-separated list of regions
    -R, --regions-file FILE           Restrict to regions listed in a file
        --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --samples LIST                List of samples to include
    -S, --samples-file FILE           File of samples to include
    -t, --targets REGION              Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE           Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
    -u, --allow-undef-tags            Print \".\" for undefined tags
    -v, --vcf-list FILE               Process multiple VCFs listed in the file
        --threads INT                 Use multithreading with INT worker threads [0]
        --verbosity INT               Verbosity level

Examples:
\tvtools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n' file.vcf.gz
\t# For more examples see http://samtools.github.io/bcftools/bcftools.html#query
";

/// Prints the usage and returns bcftools' exit status for it.
fn usage() -> ExitCode {
    eprintln!("{USAGE}");
    ExitCode::from(1)
}

fn parse_overlap(s: Option<&str>, default: Overlap, what: &str) -> anyhow::Result<Overlap> {
    match s {
        None => Ok(default),
        Some(s) => Overlap::parse(s).with_context(|| format!("Could not parse: --{what} {s}")),
    }
}

/// The record-selection options shared by every input.
struct Selection<'a> {
    threads: usize,
    regions: Option<(&'a str, ListSource)>,
    regions_overlap: Overlap,
    targets: Option<(&'a str, ListSource)>,
    targets_overlap: Overlap,
}

impl Selection<'_> {
    fn options(&self) -> anyhow::Result<SyncedOptions> {
        let mut regions = None;
        if let Some((spec, source)) = self.regions {
            let mut r = Regions::new(spec.as_bytes(), source)
                .map_err(|_| anyhow::anyhow!("Failed to read the regions: {spec}"))?;
            r.overlap = self.regions_overlap;
            regions = Some(r);
        }
        let mut targets = None;
        let mut targets_exclude = false;
        if let Some((spec, source)) = self.targets {
            let body = match spec.strip_prefix('^') {
                Some(rest) => {
                    targets_exclude = true;
                    rest
                }
                None => spec,
            };
            let mut t = Regions::new(body.as_bytes(), source)
                .map_err(|_| anyhow::anyhow!("Failed to read the targets: {spec}"))?;
            t.overlap = self.targets_overlap;
            targets = Some(t);
        }
        Ok(SyncedOptions {
            regions,
            targets,
            targets_exclude,
            apply_filters: None,
            threads: self.threads,
        })
    }
}

/// What `-s/-S` does with samples missing from the header (`--force-samples`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum UnknownSamples {
    /// Fail.
    Error,
    /// Skip them.
    Ignore,
}

fn display_name(fname: &str) -> &str {
    if fname == "-" {
        "standard input"
    } else {
        fname
    }
}

pub fn run(mut args: Args, matches: &ArgMatches) -> anyhow::Result<ExitCode> {
    if args.help {
        return Ok(usage());
    }
    if let Some(v) = &args.verbosity {
        apply_verbosity(v)?;
    }
    if args.collapse.is_some() {
        bail!("The --collapse option is obsolete, pipe through `vtools norm -c` instead.");
    }
    let filter_expr = FilterLogic::single(&args.include, &args.exclude)?;
    let format_str = format_string(&args, matches);

    let has_positional = !args.inputs.is_empty();
    let mut inputs = std::mem::take(&mut args.inputs);
    let fname = if inputs.is_empty() {
        (!io::stdin().is_terminal()).then(|| "-".to_owned())
    } else {
        Some(inputs.remove(0))
    };
    let samples = later(
        matches,
        ("samples", args.samples.as_deref()),
        ("samples_file", args.samples_file.as_deref()),
    );
    let unknown_samples = if args.force_samples {
        UnknownSamples::Ignore
    } else {
        UnknownSamples::Error
    };

    if args.list_samples {
        let Some(fname) = fname else {
            bail!("Missing the VCF file name");
        };
        let reader = SyncedReader::open(&fname, SyncedOptions::default())
            .map_err(|e| anyhow::anyhow!("Failed to read from {}: {e}", display_name(&fname)))?;
        list_columns(reader.header(), samples, unknown_samples)?;
        return Ok(ExitCode::SUCCESS);
    }

    let Some(format_str) = format_str else {
        if std::env::args().len() == 2 && fname.is_none() {
            return Ok(usage());
        }
        bail!("Error: Missing the --format option");
    };

    let out_name = args
        .output
        .as_ref()
        .map_or_else(|| "standard output".to_owned(), |p| p.display().to_string());
    let sink: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(
            File::create(path)
                .map_err(|e| anyhow::anyhow!("{}: {}", path.display(), strerror(&e)))?,
        ),
        None => Box::new(crate::output::stdout()),
    };
    let mut out = BufWriter::with_capacity(OUTPUT_BUFFER, sink);

    let selection = selection(&args, matches)?;
    let query = Query {
        format_str: &format_str,
        filter_expr,
        samples,
        unknown_samples,
        allow_undef_tags: args.allow_undef_tags,
        force_newline: !args.no_newline,
        print_filtered: args.print_filtered.as_deref(),
        out_name: &out_name,
    };

    match &args.vcf_list {
        None => {
            let Some(fname) = fname else {
                return Ok(usage());
            };
            let mut reader = open_inputs(&fname, &inputs, &selection, Message::Read)?;
            let mut prepared = query.init(&mut reader)?;
            query.print(&mut reader, &mut prepared, args.print_header, &mut out)?;
        }
        Some(list) => {
            // The positional files are read alongside each listed one.
            let mut others: Vec<String> = Vec::new();
            others.extend(fname.filter(|_| has_positional));
            others.extend(inputs);
            let listed = Listed {
                list,
                others: &others,
                selection: &selection,
                print_header: args.print_header,
            };
            query.print_listed(&listed, &mut out)?;
        }
    }
    out.flush().map_err(|e| {
        anyhow::Error::new(e).context(format!("[query_vcf] Error: cannot write to {out_name}"))
    })?;
    Ok(ExitCode::SUCCESS)
}

/// The capacity of the output buffer.
const OUTPUT_BUFFER: usize = 1 << 17;

/// The format string: `-f`, or the one `-a` makes, whichever comes later.
fn format_string(args: &Args, matches: &ArgMatches) -> Option<String> {
    let annots_format = args.annots.as_ref().map(|list| {
        let mut f = String::from("%CHROM\t%POS\t%MASK\t%REF\t%ALT\t%");
        f.push_str(&list.replace(',', "\t%"));
        f.push('\n');
        f
    });
    let last = |id: &str| last_index(matches, id);
    match (&args.format, annots_format) {
        (Some(f), Some(a)) => Some(if last("format") > last("annots") {
            f.clone()
        } else {
            a
        }),
        (Some(f), None) => Some(f.clone()),
        (None, a) => a,
    }
}

/// The regions and targets options.
fn selection<'a>(args: &'a Args, matches: &ArgMatches) -> anyhow::Result<Selection<'a>> {
    Ok(Selection {
        threads: crate::output::threads(args.threads.as_deref()).unwrap_or(0),
        regions: later(
            matches,
            ("regions", args.regions.as_deref()),
            ("regions_file", args.regions_file.as_deref()),
        ),
        regions_overlap: parse_overlap(
            args.regions_overlap.as_deref(),
            Overlap::Record,
            "regions-overlap",
        )?,
        targets: later(
            matches,
            ("targets", args.targets.as_deref()),
            ("targets_file", args.targets_file.as_deref()),
        ),
        targets_overlap: parse_overlap(
            args.targets_overlap.as_deref(),
            Overlap::Pos,
            "targets-overlap",
        )?,
    })
}

/// The inputs of `-v/--vcf-list`.
struct Listed<'a> {
    /// The file listing the inputs.
    list: &'a str,
    /// The positional files, read alongside each listed one.
    others: &'a [String],
    selection: &'a Selection<'a>,
    print_header: u8,
}

/// How a failure to open an input is worded.
#[derive(Clone, Copy)]
enum Message {
    /// `Failed to read from FILE: ...`.
    Read,
    /// `Failed to open FILE: ...` (with `-v`).
    Open,
}

/// Opens the inputs; with several, all must be indexed and are read in step, and only the
/// first one's records are printed.
fn open_inputs(
    fname: &str,
    others: &[String],
    selection: &Selection<'_>,
    message: Message,
) -> anyhow::Result<Inputs> {
    let fail = |name: &str, e: &dyn std::fmt::Display| -> anyhow::Error {
        match message {
            Message::Read => {
                anyhow::anyhow!("Failed to read from {}: {e}", display_name(name))
            }
            Message::Open => anyhow::anyhow!("Failed to open {name}: {e}"),
        }
    };
    if others.is_empty() {
        let reader =
            SyncedReader::open(fname, selection.options()?).map_err(|e| fail(fname, &e))?;
        return Ok(Inputs::single(reader));
    }
    // Several inputs need an index (HTSlib's `require_index`): open them with regions, the
    // explicit ones or, until all sequence names are known, none.
    let mut readers = Vec::with_capacity(others.len() + 1);
    for f in std::iter::once(fname).chain(others.iter().map(String::as_str)) {
        let mut options = selection.options()?;
        options.regions.get_or_insert_with(Regions::default);
        readers.push(SyncedReader::open(f, options).map_err(|e| fail(f, &e))?);
    }
    let regions = selection.options()?.regions;
    Ok(Inputs::several(readers, regions.as_ref()))
}

/// The per-input settings of a query.
struct Query<'a> {
    format_str: &'a str,
    filter_expr: Option<(&'a str, FilterLogic)>,
    samples: Option<(&'a str, ListSource)>,
    unknown_samples: UnknownSamples,
    allow_undef_tags: bool,
    force_newline: bool,
    print_filtered: Option<&'a str>,
    out_name: &'a str,
}

/// A query set up for one set of inputs.
struct Prepared {
    conv: Convert,
    filter: Option<(Box<dyn RecordFilter + Send>, FilterLogic)>,
    filter_unpack: Unpack,
}

impl Query<'_> {
    /// Sets up the samples, format, and filter for one set of inputs (bcftools' `init_data`).
    fn init(&self, reader: &mut Inputs) -> anyhow::Result<Prepared> {
        let mut order: Option<Vec<usize>> = None;
        if let Some((list, source)) = self.samples
            && list != "-"
        {
            for header in reader.headers_mut() {
                let ret = header
                    .set_samples(Some(list.as_bytes()), source)
                    .map_err(|_| anyhow::anyhow!("Error parsing the sample list"))?;
                if ret > 0 && self.unknown_samples == UnknownSamples::Error {
                    bail!(
                        "Error: sample #{ret} not found in the header, user --force-samples to proceed anyway"
                    );
                }
            }
            order = Some(sample_order(reader.header(), list, source)?);
        }
        let header = reader.header();
        let mut conv = Convert::new(header, order.as_deref(), self.format_str)?;
        if self.force_newline {
            conv.force_newline(header)?;
        }
        if self.allow_undef_tags {
            conv.set_allow_undef_tags();
        }
        let mut max_unpack = conv.max_unpack();
        let filter = match self.filter_expr {
            Some((expr, logic)) => {
                let f = filter::init(header, expr)?;
                max_unpack |= f.max_unpack();
                Some((f, logic))
            }
            None => None,
        };
        let filter_unpack = filter
            .as_ref()
            .map_or(Unpack::empty(), |(f, _)| f.max_unpack());
        if filter.is_none() || self.print_filtered.is_some() || !filter_unpack.contains(Unpack::FMT)
        {
            conv.set_header_samples();
        }
        if let Some(text) = self.print_filtered {
            conv.set_print_filtered(text.as_bytes());
        }
        reader.set_max_unpack(max_unpack);
        Ok(Prepared {
            conv,
            filter,
            filter_unpack,
        })
    }

    /// Prints the records of each file of `-v/--vcf-list`, which must all have the same
    /// samples; only the first prints the header line.
    fn print_listed(&self, listed: &Listed<'_>, out: &mut impl Write) -> anyhow::Result<()> {
        let list = listed.list;
        let Ok(fnames) = readlist::read_list(list.as_bytes(), ListSource::File) else {
            bail!("Error: failed to read {list}");
        };
        if fnames.is_empty() {
            bail!("No files in {list}?");
        }
        let mut prev_samples: Option<Vec<Box<[u8]>>> = None;
        for (i, f) in fnames.iter().enumerate() {
            let f = String::from_utf8_lossy(f).into_owned();
            let mut reader = open_inputs(&f, listed.others, listed.selection, Message::Open)?;
            let mut prepared = self.init(&mut reader)?;
            let print_header = if i == 0 { listed.print_header } else { 0 };
            // The samples are compared after -s/-S has subset them.
            if let Some(prev) = &prev_samples {
                if prev.as_slice() != reader.header().samples() {
                    bail!(
                        "Different samples in {} and {f}",
                        String::from_utf8_lossy(&fnames[i - 1])
                    );
                }
            } else {
                prev_samples = Some(reader.header().samples().to_vec());
            }
            self.print(&mut reader, &mut prepared, print_header, out)?;
        }
        Ok(())
    }

    /// Prints the header line and the records (bcftools' `query_vcf`).
    fn print(
        &self,
        reader: &mut Inputs,
        prepared: &mut Prepared,
        print_header: u8,
        out: &mut impl Write,
    ) -> anyhow::Result<()> {
        let Prepared {
            conv,
            filter,
            filter_unpack,
        } = prepared;
        let mut str = Vec::new();
        if print_header > 0 {
            if print_header > 1 {
                conv.set_no_hdr_indices();
            }
            conv.header_line(reader.header(), &mut str);
            self.write(out, &str)?;
        }
        let scope =
            if conv.max_unpack().contains(Unpack::FMT) || filter_unpack.contains(Unpack::FMT) {
                Scope::Samples
            } else {
                Scope::Site
            };
        let mut rec = Box::new(Record::new());
        let mut smpl_pass: Vec<bool> = Vec::new();
        loop {
            match reader.next_record_boxed(&mut rec) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => bail!("Error: {e}"),
            }
            // The filter and the format each decode what they read.
            let header = reader.header();
            let mut have_mask = false;
            if let Some((flt, logic)) = filter {
                let pass = flt.pass(header, &mut rec)?;
                smpl_pass.clear();
                if let Some(mask) = flt.sample_mask() {
                    smpl_pass.extend_from_slice(mask);
                    have_mask = true;
                }
                let n_sample = rec.sample_count().min(smpl_pass.len());
                let mask = have_mask.then_some(&mut smpl_pass[..n_sample]);
                if !keep(pass, *logic, scope, mask) {
                    continue;
                }
            }
            conv.line(
                header,
                &mut rec,
                have_mask.then_some(&smpl_pass[..]),
                &mut str,
            )?;
            if !str.is_empty() {
                self.write(out, &str)?;
            }
        }
        if let Some(e) = reader.take_pending_error() {
            bail!("Error: {e}");
        }
        Ok(())
    }

    fn write(&self, out: &mut impl Write, data: &[u8]) -> anyhow::Result<()> {
        out.write_all(data).map_err(|e| {
            anyhow::Error::new(e).context(format!(
                "[query_vcf] Error: cannot write to {}",
                self.out_name
            ))
        })
    }
}

/// Whether the format or the filter reads FORMAT fields.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Scope {
    /// Neither does: a record is printed or not as a whole.
    Site,
    /// Samples can be printed apart from their site.
    Samples,
}

/// Decides whether a record is printed after the filter (bcftools' `query_vcf`): `-i` can
/// include a site but exclude samples, and `-e` exclude a site but include samples.
/// `smpl_pass` holds the filter's per-sample results, when it has them.
fn keep(pass: bool, logic: FilterLogic, scope: Scope, smpl_pass: Option<&mut [bool]>) -> bool {
    match logic {
        FilterLogic::Include => pass,
        FilterLogic::Exclude => {
            if pass {
                let Some(smpl_pass) = smpl_pass else {
                    return false;
                };
                if scope == Scope::Site {
                    return false;
                }
                let mut any = false;
                for p in smpl_pass {
                    *p = !*p;
                    any |= *p;
                }
                any
            } else {
                for p in smpl_pass.into_iter().flatten() {
                    *p = true;
                }
                true
            }
        }
    }
}

/// The printed samples, in list order (bcftools' `smpl_ilist_init` with `SMPL_REORDER`;
/// header order for an exclusion list). Names are cut at the first space, as bcftools does.
fn sample_order(header: &Header, list: &str, source: ListSource) -> anyhow::Result<Vec<usize>> {
    let (negate, body) = match list.strip_prefix('^') {
        Some(rest) => (true, rest),
        None => (false, list),
    };
    if negate {
        return Ok((0..header.sample_count()).collect());
    }
    let Ok(names) = readlist::read_list(body.as_bytes(), source) else {
        bail!("Could not parse {list}");
    };
    Ok(names
        .iter()
        .filter_map(|name| {
            let end = name
                .iter()
                .enumerate()
                .position(|(i, &c)| c.is_ascii_whitespace() && !is_escaped(name, i))
                .unwrap_or(name.len());
            header.sample_id(&name[..end]).map(|i| i as usize)
        })
        .collect())
}

/// Whether the character at `i` follows an odd number of backslashes.
fn is_escaped(s: &[u8], i: usize) -> bool {
    s[..i].iter().rev().take_while(|&&c| c == b'\\').count() % 2 == 1
}

/// Prints the sample names, restricted by -s/-S (bcftools' `list_columns`).
fn list_columns(
    header: &Header,
    samples: Option<(&str, ListSource)>,
    unknown: UnknownSamples,
) -> anyhow::Result<()> {
    let mut selected: Option<(bool, std::collections::HashSet<Vec<u8>>)> = None;
    if let Some((list, source)) = samples {
        let (negate, body) = match list.strip_prefix('^') {
            Some(rest) => (true, rest),
            None => (false, list),
        };
        let Ok(names) = readlist::read_list(body.as_bytes(), source) else {
            bail!("Error: failed to read {body}");
        };
        for (i, name) in names.iter().enumerate() {
            if header.sample_id(name).is_none() && unknown == UnknownSamples::Error {
                bail!(
                    "Error: sample #{} not found in the header, user --force-samples to proceed anyway",
                    i + 1
                );
            }
        }
        selected = Some((negate, names.into_iter().collect()));
    }
    let mut out = BufWriter::new(crate::output::stdout());
    for name in header.samples() {
        let skip = match &selected {
            Some((true, set)) => set.contains(&name[..]),
            Some((false, set)) => !set.contains(&name[..]),
            None => false,
        };
        if !skip {
            out.write_all(name)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exclude_flips_the_sample_results() {
        // -e: a site that matches is printed for the samples that do not.
        let mut pass = vec![true, false, true];
        assert!(keep(
            true,
            FilterLogic::Exclude,
            Scope::Samples,
            Some(&mut pass)
        ));
        assert_eq!(pass, [false, true, false]);
        // Without per-sample results or FORMAT fields, a matching site is dropped.
        assert!(!keep(true, FilterLogic::Exclude, Scope::Samples, None));
        assert!(!keep(
            true,
            FilterLogic::Exclude,
            Scope::Site,
            Some(&mut pass)
        ));
        // A site that does not match is printed for every sample.
        assert!(keep(
            false,
            FilterLogic::Exclude,
            Scope::Samples,
            Some(&mut pass)
        ));
        assert_eq!(pass, [true, true, true]);
        assert!(!keep(
            false,
            FilterLogic::Include,
            Scope::Samples,
            Some(&mut pass)
        ));
    }

    #[test]
    fn escaped_spaces_do_not_end_sample_names() {
        assert!(is_escaped(b"a\\ b", 2));
        assert!(!is_escaped(b"a\\\\ b", 3));
        assert!(!is_escaped(b"a b", 1));
    }
}
