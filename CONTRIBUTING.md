# Contributing

Bug reports and pull requests are welcome.

## Reporting a bug

Open an issue with the version (`ntools --version` or `vtools --version`), the command line, what
you expected and why, and what happened. A small input that shows the problem helps most.
Security problems go through [SECURITY.md](SECURITY.md) instead.

## What counts as a bug

Wrong output is a bug: a record, field, count, or exit status that is not right by the format
specifications (SAM, BAM, CRAM, VCF, BCF, BGZF, and the index formats) or by what the command
documents, as well as a crash, hang, or memory blow-up. ntools and vtools follow the commands
and options of samtools, bcftools, and HTSlib 1.24, but not their mistakes: matching samtools
or bcftools does not make an output right, and differing from them does not make it wrong. A
missing or differently behaving option is worth an issue too. The [README](README.md) says
what is supported and what is left out.

## Pull requests

Before opening one, run what CI runs:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
```

The minimum supported Rust is 1.91.0, so code must build and pass clippy with it too.

A change in behaviour should come with a test that states the right output, and with a line in
[CHANGELOG.md](CHANGELOG.md) under `[Unreleased]`.

## Releases

A release tags one commit of `main` as `vX.Y.Z`; the three crates share the version. Set
`version` in `pasta/Cargo.toml`, `ntools/Cargo.toml`, and `vtools/Cargo.toml`, run `cargo build`
to update `Cargo.lock`, and rename the CHANGELOG's `[Unreleased]` to `[X.Y.Z] - YYYY-MM-DD`.
Once CI is green on that commit, push an annotated tag (`git tag -a vX.Y.Z -m "pasta X.Y.Z"`).
The release workflow checks that the tag matches the crates' version, runs the tests, and
publishes the Linux x86_64 and macOS aarch64 archives with their SHA-256 checksums.

By contributing, you agree that your contribution is licensed under the MIT License, as the rest
of the project is (see [LICENSE](LICENSE)).
