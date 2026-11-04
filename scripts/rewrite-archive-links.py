#!/usr/bin/env python3
"""Rewrites the relative links of the documents packaged in a release archive.

Usage: rewrite-archive-links.py STAGE REF REPO_PATH=ARCHIVE_PATH...

Run from the repository root. Each REPO_PATH=ARCHIVE_PATH pair says where a repository file (or,
for a directory such as `ntools`, the document that stands for it) is in the archive under
STAGE. In every packaged Markdown document, a relative link to a packaged path is rewritten to
that path relative to the document's place in the archive, and a link to anything else becomes a
link to it in the repository at REF (a tag or commit). URLs and in-page anchors are left alone.
"""

import os
import posixpath
import re
import sys

REPOSITORY = "https://github.com/DylanCosto/pasta"
LINK = re.compile(r"\]\(([^)\s]+)\)")


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)

    stage, ref = sys.argv[1], sys.argv[2]
    packaged = dict(pair.split("=", 1) for pair in sys.argv[3:])

    for repo_doc, archive_doc in packaged.items():
        if not archive_doc.endswith(".md") or not repo_doc.endswith(".md"):
            continue

        def rewrite(match):
            target = match.group(1)
            if re.match(r"[a-z]+:|#", target):
                return match.group(0)
            path, _, anchor = target.partition("#")
            suffix = f"#{anchor}" if anchor else ""
            repo_target = posixpath.normpath(posixpath.join(posixpath.dirname(repo_doc), path))

            if repo_target in packaged:
                new = posixpath.relpath(packaged[repo_target], posixpath.dirname(archive_doc) or ".")
                return f"]({new}{suffix})"

            kind = "tree" if os.path.isdir(repo_target) else "blob"
            return f"]({REPOSITORY}/{kind}/{ref}/{repo_target}{suffix})"

        file = os.path.join(stage, archive_doc)
        with open(file, encoding="utf-8") as f:
            text = f.read()
        with open(file, "w", encoding="utf-8") as f:
            f.write(LINK.sub(rewrite, text))


if __name__ == "__main__":
    main()
