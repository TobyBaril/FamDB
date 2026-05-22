#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Patch the version/date metadata of a pre-built FamDB export without
    re-running the full 6-hour export pipeline.

    Usage:
        set_ver_date.py [options] <input>

    Positional argument:
        input               Path to a FamDB directory (default) or a single
                            .h5 file when -t f is given.

    Options:
        --db-version VER    New database version string (e.g. "3.9")
        --db-date DATE      New database date in YYYY-MM-DD format
        --db-name NAME      New database name
        --db-description D  New database description (update all files)
        --file-info dump    Dump the file_info JSON from the root file to
                            <db_name>_file_info.json for manual editing
        --file-info load    Load back a previously dumped (and edited)
                            <db_name>_file_info.json into every .h5 file
        -t, --input-type    'd'/'directory' (default) or 'f'/'file'
        -n, --dry-run       Print what would change without writing anything
        -l, --log-level     Logging verbosity (default: INFO)

    What gets patched
    -----------------
    When --db-version and/or --db-date are given every .h5 file in the
    export receives:

        attrs[db_version]   <- new version
        attrs[db_date]      <- new date
        attrs[created]      <- current timestamp
        attrs[db_copyright] <- regenerated copyright block (when date changes)
        attrs[file_info]    <- embedded JSON blob's "meta.db_version" and
                               "meta.db_date" updated to match

    The file_info JSON update is essential: FamDB.__init__ reads db_version
    and db_date from that blob and cross-validates every file against it.
    Patching only the top-level attrs without updating file_info causes the
    cross-validation to fail.

SEE ALSO:
    famdb.py
    Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@systemsbiology.org>
    Anthony Gray <agray@systemsbiology.org>

LICENSE:
    This code may be used in accordance with the Creative Commons
    Zero ("CC0") public domain dedication:
    https://creativecommons.org/publicdomain/zero/1.0/

DISCLAIMER:
    This software is provided ``AS IS'' and any express or implied
    warranties, including, but not limited to, the implied warranties of
    merchantability and fitness for a particular purpose, are disclaimed.
    In no event shall the authors or the Dfam consortium members be
    liable for any direct, indirect, incidental, special, exemplary, or
    consequential damages (including, but not limited to, procurement of
    substitute goods or services; loss of use, data, or profits; or
    business interruption) however caused and on any theory of liability,
    whether in contract, strict liability, or tort (including negligence
    or otherwise) arising in any way out of the use of this software, even
    if advised of the possibility of such damage.
"""

import sys
import argparse
import datetime
import logging
import os
import re
import h5py
import json

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from famdb_globals import (
    FAMDB_ROOT_FILE_RE,
    FAMDB_COMPONENT_FILE_RE,
    META_DB_DESCRIPTION,
    META_DB_NAME,
    META_DB_COPYRIGHT,
    META_FILE_INFO,
    META_DB_VERSION,
    META_DB_DATE,
    META_CREATED,
    META_META,
    COPYRIGHT_TEXT,
)

LOGGER = logging.getLogger(__name__)


def _print_current(file_path, h5f):
    db_version = h5f.attrs.get(META_DB_VERSION, "<missing>")
    db_date = h5f.attrs.get(META_DB_DATE, "<missing>")
    meta_created = h5f.attrs.get(META_CREATED, "<missing>")
    db_name = h5f.attrs.get(META_DB_NAME, "<missing>")
    db_description = h5f.attrs.get(META_DB_DESCRIPTION, "<missing>")
    db_copyright = h5f.attrs.get(META_DB_COPYRIGHT, "<missing>")
    print(f"{file_path}:")
    print(f"  db_version  : {db_version}")
    print(f"  db_date     : {db_date}")
    print(f"  created     : {meta_created}")
    print(f"  db_name     : {db_name}")
    print(f"  description : {db_description[:60]}{'...' if len(str(db_description)) > 60 else ''}")
    print(f"  copyright   : {str(db_copyright)[:60]}...")


def update_file(file_path, new_creation_time, args, dry_run=False):
    """Open one .h5 file and apply the requested metadata patches."""
    open_mode = "r" if dry_run else "r+"
    with h5py.File(file_path, mode=open_mode) as h5f:
        _print_current(file_path, h5f)

        if dry_run:
            _print_pending(args)
            return

        # --- db_version ---
        if args.db_version:
            h5f.attrs[META_DB_VERSION] = args.db_version
            print(f"    ** new db_version  : {args.db_version}")

        # --- db_name ---
        if args.db_name:
            h5f.attrs[META_DB_NAME] = args.db_name
            print(f"    ** new db_name     : {args.db_name}")

        # --- db_description ---
        if args.db_description:
            h5f.attrs[META_DB_DESCRIPTION] = args.db_description
            print(f"    ** new description : {args.db_description[:60]}")

        # --- db_date (also regenerates copyright and bumps created) ---
        if args.db_date:
            year_match = re.match(r"^(\d{4})-\d{2}-\d{2}$", args.db_date)
            if not year_match:
                raise ValueError(
                    f"Date must be YYYY-MM-DD, got: {args.db_date!r}"
                )
            db_year = year_match.group(1)
            effective_version = args.db_version or h5f.attrs.get(META_DB_VERSION, "")
            copyright_text = COPYRIGHT_TEXT % (
                db_year,
                effective_version,
                args.db_date,
            )
            h5f.attrs[META_DB_DATE] = args.db_date
            h5f.attrs[META_CREATED] = new_creation_time
            h5f.attrs[META_DB_COPYRIGHT] = copyright_text
            print(f"    ** new db_date     : {args.db_date}")
            print(f"    ** new created     : {new_creation_time}")
            print(f"    ** new copyright   : {copyright_text[:60]}...")

        # --- file_info JSON: update meta.db_version / meta.db_date ---
        # This blob is read by FamDB.__init__ and cross-validated across
        # all component files.  If only the top-level attrs are patched the
        # load will fail with a "Files From Different Partitioning Runs" error.
        if args.db_version or args.db_date:
            raw = h5f.attrs.get(META_FILE_INFO)
            if raw is not None:
                info_obj = json.loads(raw)
                meta_section = info_obj.get(META_META, {})
                if args.db_version:
                    meta_section[META_DB_VERSION] = args.db_version
                if args.db_date:
                    meta_section[META_DB_DATE] = args.db_date
                info_obj[META_META] = meta_section
                h5f.attrs[META_FILE_INFO] = json.dumps(info_obj)
                print(f"    ** file_info JSON  : meta.db_version/db_date updated")

        # --- dump file_info JSON to a file ---
        db_name = h5f.attrs.get(META_DB_NAME, "famdb")
        dump_name = f"{db_name}_file_info.json"
        if args.file_info == "dump":
            raw = h5f.attrs.get(META_FILE_INFO)
            if raw is None:
                LOGGER.warning(f"No {META_FILE_INFO} attribute found in {file_path}")
            else:
                info_obj = json.loads(raw)
                with open(dump_name, "w") as outfile:
                    json.dump(info_obj, outfile, indent=4)
                print(f"    ** file_info dumped: {dump_name}")

        # --- load file_info JSON from a file ---
        if args.file_info == "load":
            if not os.path.exists(dump_name):
                raise FileNotFoundError(
                    f"Cannot load file_info: {dump_name!r} not found. "
                    f"Run with --file-info dump first."
                )
            with open(dump_name, "r") as infile:
                new_info = json.load(infile)
            h5f.attrs[META_FILE_INFO] = json.dumps(new_info)
            print(f"    ** file_info loaded: {dump_name}")


def _print_pending(args):
    """Print what would be changed in dry-run mode."""
    changes = []
    if args.db_version:
        changes.append(f"  would set db_version  -> {args.db_version}")
    if args.db_date:
        changes.append(f"  would set db_date     -> {args.db_date}")
        changes.append(f"  would set created     -> <now>")
        changes.append(f"  would regenerate copyright")
    if args.db_name:
        changes.append(f"  would set db_name     -> {args.db_name}")
    if args.db_description:
        changes.append(f"  would set description -> {args.db_description[:60]}")
    if args.db_version or args.db_date:
        changes.append(f"  would update file_info JSON meta section")
    if args.file_info:
        changes.append(f"  would {args.file_info} file_info JSON")
    if changes:
        print("  [DRY RUN] pending changes:")
        for c in changes:
            print(c)
    else:
        print("  [DRY RUN] no changes requested")


def collect_h5_files(input_path, input_type):
    """
    Return a list of .h5 file paths to process.

    For directory mode: discovers the root file (*.0.h5) plus all component
    files matching the v3 naming convention.  Non-famdb .h5 files are skipped.
    For file mode: returns [input_path] if it looks like a famdb file.
    """
    if input_type in ("f", "file"):
        basename = os.path.basename(input_path)
        if FAMDB_ROOT_FILE_RE.match(basename) or FAMDB_COMPONENT_FILE_RE.match(basename):
            return [input_path]
        # Also accept the old single-partition pattern *.N.h5
        if re.match(r"\S+\.\d+\.h5$", basename):
            return [input_path]
        LOGGER.warning(
            f"File {input_path!r} does not match a known FamDB filename pattern — "
            "processing anyway."
        )
        return [input_path]

    # directory mode
    if not os.path.isdir(input_path):
        raise NotADirectoryError(f"Expected a directory, got: {input_path!r}")

    files = []
    for filename in sorted(os.listdir(input_path)):
        if FAMDB_ROOT_FILE_RE.match(filename) or FAMDB_COMPONENT_FILE_RE.match(filename):
            files.append(os.path.join(input_path, filename))
    if not files:
        LOGGER.warning(f"No FamDB .h5 files found in {input_path!r}")
    return files


def main():
    logging.basicConfig(format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(
        description="Patch version/date metadata in a FamDB export."
    )
    parser.add_argument("-l", "--log-level", default="INFO")
    parser.add_argument("--db-version", metavar="VER",
                        help="New database version string (e.g. '3.9')")
    parser.add_argument("--db-date", metavar="YYYY-MM-DD",
                        help="New database date")
    parser.add_argument("--db-name", metavar="NAME",
                        help="New database name")
    parser.add_argument("--db-description", metavar="DESC",
                        help="New database description")
    parser.add_argument("--file-info", choices=("dump", "load"),
                        help="Dump or load the file_info JSON for manual editing")
    parser.add_argument("-t", "--input-type",
                        choices=("f", "file", "d", "directory"),
                        default="d",
                        help="Input type: 'd'/'directory' (default) or 'f'/'file'")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="Show what would change without writing")
    parser.add_argument("input",
                        help="Path to FamDB directory or single .h5 file")

    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    if not any([args.db_version, args.db_date, args.db_name,
                args.db_description, args.file_info]):
        parser.error(
            "Nothing to do — specify at least one of: "
            "--db-version, --db-date, --db-name, --db-description, --file-info"
        )

    new_creation_time = str(datetime.datetime.now())

    try:
        h5_files = collect_h5_files(args.input, args.input_type)
    except (NotADirectoryError, FileNotFoundError) as exc:
        LOGGER.error(str(exc))
        sys.exit(1)

    if not h5_files:
        LOGGER.error("No .h5 files to process.")
        sys.exit(1)

    errors = []
    for file_path in h5_files:
        try:
            update_file(file_path, new_creation_time, args, dry_run=args.dry_run)
        except Exception as exc:
            LOGGER.error(f"Failed to update {file_path}: {exc}")
            errors.append(file_path)

    if errors:
        LOGGER.error(f"{len(errors)} file(s) could not be updated: {errors}")
        sys.exit(1)

    if not args.dry_run:
        print(f"\nDone. Updated {len(h5_files)} file(s).")


if __name__ == "__main__":
    main()
