#!/usr/bin/env python3
"""Writes the licence texts of the crates linked into ntools and vtools.

Usage: third-party-licenses.py OUTPUT

Reads `cargo metadata` for the workspace (Linux x86_64, default features), follows the normal
(non-dev, non-build) dependencies of ntools and vtools, and writes, for each crate, its name,
version and licence expression, and every licence or notice file it ships: LICENSE*, LICENCE*,
COPYING*, NOTICE*, UNLICENSE, and the same names one directory down, which is where `-sys`
crates keep the C library they bundle (libbz2-rs-sys's, say). A text identical to one already
written (apart from whitespace) is written once and referred to after that. Run it after
`cargo build`, so the crates' sources are in Cargo's registry. Exits nonzero if a crate ships
no licence file.
"""

import json
import pathlib
import subprocess
import sys

NAMES = ("license", "licence", "copying", "notice", "unlicense")


def licence_files(directory):
    found = []
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.name.lower().startswith(NAMES):
            found.append(path)
        elif path.is_dir() and not path.name.startswith("."):
            found.extend(
                p
                for p in sorted(path.iterdir())
                if p.is_file() and p.name.lower().startswith(NAMES)
            )
    return found


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)

    metadata = json.loads(
        subprocess.run(
            [
                "cargo", "metadata", "--format-version", "1", "--locked",
                "--filter-platform", "x86_64-unknown-linux-gnu",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )

    packages = {p["id"]: p for p in metadata["packages"]}
    nodes = {n["id"]: n for n in metadata["resolve"]["nodes"]}
    members = set(metadata["workspace_members"])

    stack = [i for i in members if packages[i]["name"] in ("ntools", "vtools")]
    linked = set()

    while stack:
        package_id = stack.pop()
        if package_id in linked:
            continue
        linked.add(package_id)
        for dep in nodes[package_id]["deps"]:
            if any(kind["kind"] is None for kind in dep["dep_kinds"]):
                stack.append(dep["pkg"])

    crates = sorted(
        (packages[i] for i in linked - members), key=lambda p: (p["name"], p["version"])
    )

    missing = []
    # Whitespace-normalized text -> the crate and file it was first written for.
    written = {}

    with open(sys.argv[1], "w", encoding="utf-8") as out:
        out.write(
            "Third-party software in the ntools and vtools binaries\n\n"
            "These crates are compiled into the binaries. Each entry gives the crate, its\n"
            "licence expression as published, and the licence and notice files it ships.\n"
            "pasta's own licence is in LICENSE, and the projects it derives from are in NOTICE.\n"
        )

        for crate in crates:
            directory = pathlib.Path(crate["manifest_path"]).parent
            files = licence_files(directory)
            out.write("\n" + "=" * 78 + "\n")
            out.write(f"{crate['name']} {crate['version']}\n")
            out.write(f"Licence: {crate.get('license') or 'see the files below'}\n")
            if crate.get("repository"):
                out.write(f"Source: {crate['repository']}\n")

            if not files:
                missing.append(f"{crate['name']} {crate['version']}")

            for path in files:
                text = path.read_text(encoding="utf-8", errors="replace").rstrip()
                name = path.relative_to(directory)
                key = " ".join(text.split())
                out.write("\n" + "-" * 78 + "\n")
                if key in written:
                    out.write(f"{name}: the same text as {written[key]} above.\n")
                else:
                    written[key] = f"{name} of {crate['name']} {crate['version']}"
                    out.write(f"{name}:\n\n{text}\n")

    if missing:
        sys.exit("no licence file shipped by: " + ", ".join(missing))


if __name__ == "__main__":
    main()
