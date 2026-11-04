# Security

pasta, ntools, and vtools read files that may come from anywhere, so a crash, hang, or memory
blow-up on a crafted SAM, BAM, CRAM, VCF, BCF, FASTA, or index file is a bug worth reporting.

## Reporting a vulnerability

Report it privately through GitHub's private vulnerability reporting: the **Security** tab of
this repository, then **Report a vulnerability**. Please include the version (`ntools --version`
or `vtools --version`, or the commit), the command line, and an input that triggers it, made as
small as you can.

Please don't open a public issue for a vulnerability until it is fixed. Bugs that are not
security problems (wrong results, missing options) can go in public issues.

## Supported versions

Fixes go into the next release; only the latest release is supported.
