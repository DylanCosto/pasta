//! Expressions evaluated on small records, with the results bcftools 1.24 gives.

use pasta::variant::{Header, Record, vcf};

use super::Filter;

const HEADER: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FILTER=<ID=q10,Description=\"Low quality\">\n\
##contig=<ID=1>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"f\">\n\
##INFO=<ID=STR,Number=.,Type=String,Description=\"s\">\n\
##INFO=<ID=FL,Number=0,Type=Flag,Description=\"f\">\n\
##INFO=<ID=CNT,Number=A,Type=Integer,Description=\"c\">\n\
##INFO=<ID=DOT,Number=.,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
##FORMAT=<ID=XD,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

const LINES: [&str; 3] = [
    "1\t100\trs1;rs2\tA\tG\t50\tPASS\tDP=10;AF=0.25;STR=x,yy;FL\tGT:AD:XD\t0/1:5,5:10\t1/1:0,9:9\t./.:.:.",
    "1\t200\t.\tAC\tA,ACC\t.\tq10\tAF=0.1,0.5\tGT:AD\t0|0:3,0,0\t1|2:0,4,4\t0:7,.,.",
    "1\t300\t.\tG\t<DEL>\t5\t.\t.\tGT\t0/0\t0/0\t0/0",
];

/// Evaluates `expr` on each line, returning the site results and the per-sample results.
fn eval(expr: &str) -> Vec<(bool, Option<Vec<bool>>)> {
    eval_lines(expr, &LINES)
}

/// Evaluates `expr` on each of `lines`, as [`eval`] does.
fn eval_lines(expr: &str, lines: &[&str]) -> Vec<(bool, Option<Vec<bool>>)> {
    let mut header = Header::parse(HEADER.as_bytes()).expect("the header parses");
    let mut filter = Filter::init(&header, expr).expect("the expression compiles");
    let mut rec = Record::new();
    lines
        .iter()
        .map(|line| {
            let mut line = line.as_bytes().to_vec();
            vcf::parse(&mut header, &mut line, &mut rec).expect("the line parses");
            let pass = filter
                .test(&header, &mut rec)
                .expect("the expression evaluates");
            (pass, filter.sample_pass().map(<[bool]>::to_vec))
        })
        .collect()
}

fn sites(expr: &str) -> Vec<bool> {
    eval(expr).into_iter().map(|(p, _)| p).collect()
}

#[test]
fn site_values() {
    assert_eq!(sites("QUAL>10"), [true, false, false]);
    assert_eq!(sites("QUAL=\".\""), [false, true, false]);
    assert_eq!(sites("AF[1]>0.4"), [false, true, false]);
    assert_eq!(sites("ID~\"^rs2\""), [true, false, false]);
    assert_eq!(sites("ID!=\"rs1\""), [false, true, true]);
    assert_eq!(sites("FILTER=\"PASS\""), [true, false, false]);
    assert_eq!(sites("FILTER=\".\""), [false, false, true]);
    assert_eq!(sites("TYPE=\"indel\""), [false, true, false]);
    assert_eq!(sites("TYPE~\"other\""), [false, false, true]);
    assert_eq!(sites("ILEN=-1"), [false, true, false]);
    assert_eq!(sites("N_ALT=2"), [false, true, false]);
    assert_eq!(sites("STR=\"yy\""), [true, false, false]);
    assert_eq!(sites("FL=1"), [true, false, false]);
    assert_eq!(sites("ALT=\"ACC\""), [false, true, false]);
    assert_eq!(sites("strlen(REF)=2"), [false, true, false]);
}

#[test]
fn computed_values() {
    assert_eq!(sites("AN=4"), [true, false, false]);
    assert_eq!(sites("AC[*]=3"), [true, false, false]);
    assert_eq!(sites("MAF<0.3"), [true, true, true]);
    assert_eq!(sites("N_MISSING=1"), [true, false, false]);
    assert_eq!(sites("F_MISSING>0.3"), [true, false, false]);
}

#[test]
fn functions() {
    assert_eq!(sites("SUM(FMT/AD)=19"), [true, false, false]);
    assert_eq!(sites("MAX(FMT/AD)=9"), [true, false, false]);
    assert_eq!(sites("MEDIAN(FMT/AD)=4"), [false, false, false]);
    assert_eq!(sites("COUNT(GT=\"het\")=1"), [true, true, false]);
    assert_eq!(sites("N_PASS(GT=\"hom\")=2"), [false, false, false]);
    assert_eq!(sites("N_PASS(GT=\"hom\")=3"), [false, false, true]);
    // bcftools passes the third record too, from the previous record's AD values.
    assert_eq!(sites("binom(FMT/AD)=1"), [true, true, false]);
    assert_eq!(sites("ABS(-INFO/DP)=10"), [true, false, false]);
}

#[test]
fn arithmetic_is_right_associative() {
    // bcftools pops operators only of higher precedence, so 10-3-2 is 10-(3-2).
    assert_eq!(sites("10-3-2=9"), [true, true, true]);
    assert_eq!(sites("1+2*3=7"), [true, true, true]);
}

#[test]
fn per_sample_results() {
    let r = eval("GT=\"het\"");
    assert_eq!(r[0].1.as_deref(), Some(&[true, false, false][..]));
    assert_eq!(r[1].1.as_deref(), Some(&[false, true, false][..]));
    let r = eval("FMT/XD>9 && GT=\"alt\"");
    assert_eq!(r[0].1.as_deref(), Some(&[true, true, false][..]));
    let r = eval("FMT/XD>9 & GT=\"alt\"");
    assert_eq!(r[0].1.as_deref(), Some(&[true, false, false][..]));
    let r = eval("QUAL>10");
    assert_eq!(r[0].1, None);
}

/// Records with mixed ploidy and a two-value `AF`.
const PLOIDY_LINES: [&str; 2] = [
    "1\t400\t.\tA\tC,G\t9\t.\tAF=0.5,0.1\tGT:AD:XD\t0/1:5,6,0:10\t2:1,0,7:4\t1/2:3,3,3:.",
    "1\t500\t.\tA\tC\t9\t.\tAF=0.2\tGT:AD:XD\t0/1/1:5,6:1\t1:1,5:4\t0/0:3,5:7",
];

#[test]
fn per_sample_arithmetic_with_a_vector_uses_its_first_value() {
    let r = eval_lines("FMT/XD*AF>2", &PLOIDY_LINES);
    assert_eq!(r[0], (true, Some(vec![true, false, false])));
    assert_eq!(r[1], (false, Some(vec![false, false, false])));
}

#[test]
fn gt_subscripted_comparisons_are_symmetric() {
    // bcftools gives sample B of the second record for 5!=AD[*:GT] (its padding read as
    // missing), but not for AD[*:GT]!=5.
    for expr in ["FMT/AD[*:GT]!=5", "5!=FMT/AD[*:GT]"] {
        let r = eval_lines(expr, &PLOIDY_LINES);
        assert_eq!(r[1].1.as_deref(), Some(&[true, false, true][..]), "{expr}");
    }
}

#[test]
fn missing_genotypes_are_counted_per_sample() {
    // Diploid genotypes (the 8-bit fast path) and a haploid one padded with the vector end.
    let lines = [
        "1\t1\t.\tA\tC\t.\t.\t.\tGT\t./.\t0/.\t./1",
        "1\t2\t.\tA\tC\t.\t.\t.\tGT\t.|.\t0|1\t1|1",
        "1\t3\t.\tA\tC\t.\t.\t.\tGT\t.\t0\t0/1",
    ];
    let count = |n: u8| sites_of(eval_lines(&format!("N_MISSING={n}"), &lines));
    assert_eq!(count(3), [true, false, false]);
    assert_eq!(count(1), [false, true, true]);
    assert_eq!(
        sites_of(eval_lines("F_MISSING>0.5", &lines)),
        [true, false, false]
    );
}

fn sites_of(results: Vec<(bool, Option<Vec<bool>>)>) -> Vec<bool> {
    results.into_iter().map(|(p, _)| p).collect()
}

/// Records with missing integer values: alone, before a value, and after one.
const MISSING_INT_LINES: [&str; 4] = [
    "1\t1\t.\tA\tC\t.\t.\tCNT=.;DOT=.\tGT:AD\t0/1:.\t0/1:.,5\t0/1:3,.",
    "1\t2\t.\tA\tC,G\t.\t.\tCNT=.,1;DOT=.,1\tGT:AD\t0/1:.,.,.\t0/2:2,.,7\t1/1:3,4,5",
    "1\t3\t.\tA\tC,G\t.\t.\tCNT=1,.;DOT=1,.\tGT:AD\t0/1:1,2,3\t0/2:1,2,3\t1/1:1,2,3",
    "1\t4\t.\tA\tC\t.\t.\tCNT=1;DOT=1\tGT:AD\t0/1:2,5\t0/1:3,4\t0/1:1,1",
];

#[test]
fn missing_integers_are_missing_values() {
    let missing_sites = |expr: &str| sites_of(eval_lines(expr, &MISSING_INT_LINES));
    // bcftools reads these `.` as the numbers -128 or -2147483648, so they are all true there.
    for expr in [
        "CNT<0",
        "SUM(CNT)<0",
        "MIN(CNT)<0",
        "CNT[*]<0",
        "CNT[0]<0",
        "CNT[1]<0",
        "DOT<0",
        "FMT/AD<0",
        "FMT/AD[*:GT]<0",
    ] {
        assert_eq!(missing_sites(expr), [false; 4], "{expr}");
    }
    assert_eq!(missing_sites("SUM(CNT)=1"), [false, true, true, true]);
    assert_eq!(missing_sites("CNT=\".\""), [true, true, true, false]);
    assert_eq!(
        missing_sites("FMT/AD[*:GT]=\".\""),
        [true, true, false, false]
    );
    // An index past a single value reads nothing; bcftools reads the value at any index.
    assert_eq!(missing_sites("CNT[1]=1"), [false, true, false, false]);
    assert_eq!(sites("AF[1]=0.25"), [false; 3]);
}

#[test]
fn errors() {
    let header = Header::parse(HEADER.as_bytes()).expect("the header parses");
    let message = |expr: &str| {
        Filter::init(&header, expr)
            .expect_err("the expression is rejected")
            .to_string()
    };
    assert_eq!(
        message("NOSUCH>1"),
        "[filter.c:3615 filters_init1] Error: the tag \"NOSUCH\" is not defined in the VCF header"
    );
    assert_eq!(message("\"abc"), "Missing quotes in: \"abc");
    assert_eq!(
        message("(QUAL>1"),
        "Could not parse the expression: [(QUAL>1]"
    );
    assert_eq!(
        message("FMT/AD[0]>1"),
        "The FORMAT tag AD can have multiple subfields, run as AD[sample:subfield]"
    );
    assert_eq!(
        message("GT=\"het\" || GT="),
        "Comparison between samples is not supported, sorry!"
    );
    assert_eq!(
        message("FILTER~"),
        "[filter.c:4194 filter_init_] Could not parse the expression, an unquoted string value perhaps? FILTER~"
    );
    assert_eq!(
        message("TYPE=\"foo\""),
        "The type \"foo\" not recognised: TYPE=\"foo\""
    );
}
