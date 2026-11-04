#!/usr/bin/env bash
# Packages built ntools and vtools binaries as a release archive.
#
# Usage: package-release.sh VERSION [OUT_DIR]
#
# Run from the repository root after `cargo build --release --locked -p ntools -p vtools`, on
# Linux x86_64 or macOS aarch64. Writes OUT_DIR (default dist)/pasta-VERSION-PLATFORM.tar.gz
# (PLATFORM is x86_64-linux or aarch64-macos) and its .sha256. The archive holds the two
# binaries, README.md, CHANGELOG.md, LICENSE, NOTICE, the crates' READMEs as docs/ntools.md,
# docs/vtools.md and docs/pasta.md, THIRD_PARTY_LICENSES.txt (the licences of the crates linked
# in, from scripts/third-party-licenses.py), and licenses/ with the licence texts of the projects
# ntools and vtools port (samtools, HTSlib, htscodecs, bcftools with its edlib, and bcftools'
# VariantKey), taken from their 1.24 release tarballs, downloaded unless already in
# $DOWNLOAD_DIR.
# Links in the packaged documents to files left out point at the repository at VERSION when it
# is a tag (vX.Y.Z), or else at the commit being packaged.
#
# On Linux, the binaries must not need a glibc newer than $MAX_GLIBC (default 2.28, RHEL 8's), so
# that they run on RHEL 8 and later; the script fails otherwise.
set -euo pipefail

version=$1
out=${2:-dist}
max_glibc=${MAX_GLIBC:-2.28}
upstream=1.24
downloads=${DOWNLOAD_DIR:-$(mktemp -d)}
case "$(uname -s)-$(uname -m)" in
    Linux-x86_64) platform=x86_64-linux ;;
    Darwin-arm64) platform=aarch64-macos ;;
    *) echo "unsupported platform: $(uname -s) $(uname -m)" >&2; exit 1 ;;
esac
name=pasta-$version-$platform
stage=$out/$name

if [ "$platform" = x86_64-linux ]; then
    for binary in target/release/ntools target/release/vtools; do
        needed=$(objdump -T "$binary" | grep -o 'GLIBC_[0-9.]*' | sed 's/GLIBC_//' | sort -V | tail -1)
        if [ "$(printf '%s\n%s\n' "$needed" "$max_glibc" | sort -V | tail -1)" != "$max_glibc" ]; then
            echo "$binary needs glibc $needed, newer than $max_glibc" >&2
            exit 1
        fi
        echo "$binary needs glibc $needed"
    done
fi

rm -rf "$stage"
mkdir -p "$stage/licenses"
cp target/release/ntools target/release/vtools README.md CHANGELOG.md LICENSE NOTICE "$stage/"
# The crates' READMEs go under docs/, since ntools and vtools name the binaries.
mkdir -p "$stage/docs"
for crate in ntools vtools pasta; do
    cp "$crate/README.md" "$stage/docs/$crate.md"
done
"${PYTHON:-python3}" scripts/third-party-licenses.py "$stage/THIRD_PARTY_LICENSES.txt"

# The ported projects' own licence texts, which NOTICE refers to.
for tool in samtools bcftools; do
    tarball=$downloads/$tool-$upstream.tar.bz2
    if [ ! -f "$tarball" ]; then
        curl -fsSL -o "$tarball" "https://github.com/samtools/$tool/releases/download/$upstream/$tool-$upstream.tar.bz2"
    fi
    tar -xjf "$tarball" -C "$downloads" "$tool-$upstream/LICENSE"
    cp "$downloads/$tool-$upstream/LICENSE" "$stage/licenses/$tool-LICENSE"
done
bcftools_src=bcftools-$upstream
tar -xjf "$downloads/$bcftools_src.tar.bz2" -C "$downloads" "$bcftools_src/htslib-$upstream/LICENSE" \
    "$bcftools_src/htslib-$upstream/htscodecs/LICENSE.md" "$bcftools_src/variantkey.h"
cp "$downloads/$bcftools_src/htslib-$upstream/LICENSE" "$stage/licenses/htslib-LICENSE"
cp "$downloads/$bcftools_src/htslib-$upstream/htscodecs/LICENSE.md" "$stage/licenses/htscodecs-LICENSE.md"
# VariantKey's licence is the comment block at the top of its header.
sed -n '/^\/\/ LICENSE$/,/^\/\/ THE SOFTWARE\.$/s|^// \{0,1\}||p' "$downloads/$bcftools_src/variantkey.h" \
    >"$stage/licenses/variantkey-LICENSE"
grep -q "Nicola Asuni" "$stage/licenses/variantkey-LICENSE"
# vtools ports bcftools' edlib.c, whose licence NOTICE says bcftools' LICENSE holds.
grep -q "License for edlib" "$stage/licenses/bcftools-LICENSE"

if [[ $version =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    ref=$version
else
    ref=${GITHUB_SHA:-$(git rev-parse HEAD)}
fi
"${PYTHON:-python3}" scripts/rewrite-archive-links.py "$stage" "$ref" \
    README.md=README.md CHANGELOG.md=CHANGELOG.md LICENSE=LICENSE NOTICE=NOTICE \
    ntools/README.md=docs/ntools.md vtools/README.md=docs/vtools.md pasta/README.md=docs/pasta.md \
    ntools=docs/ntools.md vtools=docs/vtools.md pasta=docs/pasta.md

if command -v sha256sum >/dev/null; then
    sha256=(sha256sum)
else
    sha256=(shasum -a 256)
fi
(cd "$out" && tar -czf "$name.tar.gz" "$name" && "${sha256[@]}" "$name.tar.gz" >"$name.tar.gz.sha256")
echo "wrote $out/$name.tar.gz"
