#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge RepBase RepeatMasker Edition (RMRB) families into locally-installed
FamDB curated-consensus partitions.

RepBase is distributed as two EMBL files:
  RMRBMeta.embl  - metadata (taxonomy, classification, type/subtype)
  RMRBSeqs.embl  - sequences (must be obtained separately from GIRI/RepBase)

This script combines them into a single RMRB.embl (cached for reuse) and then
appends any families not already present into the appropriate CC partitions.

STATE TRACKING
  A JSON file (.repbase_merge_state.json) in the FamDB directory records which
  CC partition files have been merged and what version of RepBase was used.
  Re-running the script is safe: existing entries are skipped, only newly
  installed partitions (or a changed RepBase file) trigger fresh work.

INTERRUPTED MERGES
  A 'merge.working' sentinel file is written before any modifications and
  removed on success.  If it is still present on the next run, the FamDB
  files may be in an inconsistent state.  Delete all files in the FamDB
  directory, re-download the partitions, and re-run this script.

Usage:
    merge_repbase.py -i <famdb_dir>
                     [--meta RMRBMeta.embl] [--seqs RMRBSeqs.embl]
                     [--combined RMRB.embl] [--dup RMRB_DUP.txt]
                     [--name NAME] [--description DESC]
                     [--force] [-l LOG_LEVEL]

DEFAULT FILE LOCATIONS
  The script first looks for source files in the Libraries/ directory
  (sibling of the utils/ directory containing this script):
    Libraries/RMRBMeta.embl
    Libraries/RMRBSeqs.embl
    Libraries/RMRB_DUP.txt

  If a file is not found there, an error is raised and the corresponding
  command-line flag (--meta, --seqs, --dup) must be used to supply the path.
"""

import argparse
import datetime
import json
import logging
import os
import re
import sys

# Allow importing famdb_classes from the project root when run as a script
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from famdb_classes import FamDB
from famdb_globals import COMPONENT_CC, FAMDB_COMPONENT_FILE_RE

LOGGER = logging.getLogger(__name__)

STATE_FILE = ".repbase_merge_state.json"
SENTINEL_FILE = "merge.working"


# ---------------------------------------------------------------------------
# EMBL helpers
# ---------------------------------------------------------------------------

def _get_embl_version(filepath):
    """Return the RELEASE string from the header of an EMBL file, or ''."""
    try:
        with open(filepath) as fh:
            for i, line in enumerate(fh):
                m = re.match(r"^(?:CC|##)\s+RELEASE\s+(\S+);", line)
                if m:
                    return m.group(1)
                if i > 60:
                    break
    except OSError:
        pass
    return ""


def _read_embl_sequences(seqs_file):
    """
    Parse an EMBL file that contains sequence data (RMRBSeqs.embl).
    Returns {accession: sequence_string}.
    """
    sequences = {}
    current_id = None
    current_seq = None
    in_seq = False

    with open(seqs_file) as fh:
        for line in fh:
            if line.startswith("ID"):
                m = re.match(r"ID\s+(\S+)", line)
                if m:
                    current_id = m.group(1).rstrip(";")
                    current_seq = ""
                    in_seq = False
            elif line.startswith("SQ"):
                in_seq = True
            elif line.startswith("//"):
                if current_id is not None and current_seq is not None:
                    sequences[current_id] = current_seq
                current_id = None
                current_seq = None
                in_seq = False
            elif in_seq and current_seq is not None:
                current_seq += re.sub(r"[^A-Za-z]", "", line)

    return sequences


def combine_rmrb_files(meta_file, seqs_file, output_file):
    """
    Combine RMRBMeta.embl (metadata) with RMRBSeqs.embl (sequences) into a
    single RMRB.embl.  Returns the RepBase version string found in seqs_file.
    """
    LOGGER.info(f"Reading sequences from {seqs_file} ...")
    sequences = _read_embl_sequences(seqs_file)
    LOGGER.info(f"  Read {len(sequences):,} sequences")

    version = _get_embl_version(seqs_file) or _get_embl_version(meta_file)

    header = (
        "CC ****************************************************************\n"
        "CC                                                                *\n"
        "CC   RepBase RepeatMasker Edition                                 *\n"
        "CC    Please refer to GIRI (https://www.girinst.org/) for         *\n"
        "CC    detailed copyright and licensing restrictions.              *\n"
        "CC                                                                *\n"
        f"CC    RELEASE {version};{' ' * max(0, 35 - len(version))}*\n"
        "CC                                                                *\n"
        "CC   RepeatMasker software and maintenance are currently          *\n"
        "CC   funded by an NIH/NHGRI R01 grant HG02939-01 to Arian Smit.  *\n"
        "CC                                                                *\n"
        "CC ****************************************************************\n"
    )

    LOGGER.info(f"Reading metadata from {meta_file} ...")
    merged = 0
    missing = 0

    with open(output_file, "w") as out:
        out.write(header)

        current_id = None
        meta_lines = []
        in_record = False

        with open(meta_file) as mfh:
            for line in mfh:
                if line.startswith("ID"):
                    current_id = None
                    meta_lines = []
                    in_record = True
                    m = re.match(r"ID\s+(\S+)", line)
                    if m:
                        current_id = m.group(1).rstrip(";")
                    meta_lines.append(line)

                elif line.startswith("//") and in_record:
                    if current_id:
                        seq = sequences.get(current_id)
                        if seq:
                            # Fix "???" placeholder with actual length
                            meta_lines[0] = re.sub(r"\?\?\?", str(len(seq)), meta_lines[0])
                            for ml in meta_lines:
                                out.write(ml)
                            out.write("CC   Source: RepBase RepeatMasker Edition\n")
                            out.write(f"SQ   Sequence {len(seq)} BP;\n")
                            # Write in 60-char lines (no position numbers needed;
                            # read_embl_families strips all non-alpha chars)
                            for i in range(0, len(seq), 60):
                                out.write("     " + seq[i:i + 60] + "\n")
                            out.write("//\n")
                            merged += 1
                        else:
                            LOGGER.warning(f"No sequence found for {current_id}; skipping")
                            missing += 1
                    current_id = None
                    meta_lines = []
                    in_record = False

                elif in_record:
                    meta_lines.append(line)

    LOGGER.info(f"  Combined {merged:,} records ({missing} skipped - no sequence)")
    return version


# ---------------------------------------------------------------------------
# State file helpers
# ---------------------------------------------------------------------------

def _state_path(db_dir):
    return os.path.join(db_dir, STATE_FILE)


def _sentinel_path(db_dir):
    return os.path.join(db_dir, SENTINEL_FILE)


def load_state(db_dir):
    path = _state_path(db_dir)
    if os.path.exists(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except Exception as exc:
            LOGGER.warning(f"Could not read state file {path}: {exc}")
    return {"repbase_version": "", "rmrb_file": "", "rmrb_size": 0, "merged_cc": {}}


def save_state(db_dir, state):
    path = _state_path(db_dir)
    with open(path, "w") as fh:
        json.dump(state, fh, indent=2)
    LOGGER.debug(f"State written to {path}")


# ---------------------------------------------------------------------------
# Partition scanning
# ---------------------------------------------------------------------------

def get_cc_partition_files(db_dir):
    """Return {filename: size} for every CC component file in db_dir."""
    result = {}
    for fname in os.listdir(db_dir):
        m = FAMDB_COMPONENT_FILE_RE.match(fname)
        if m and m.group(2) == "curated" and m.group(3) == "consensus":
            full = os.path.join(db_dir, fname)
            result[fname] = os.path.getsize(full)
    return result


def check_needs_merge(state, current_cc, rmrb_file):
    """
    Returns (needs_merge: bool, reason: str).

    A merge is needed when:
      - any CC partition file is absent from the state (newly downloaded), OR
      - the RMRB file path/size differs from what was last used.
    """
    new_parts = [f for f in current_cc if f not in state.get("merged_cc", {})]
    if new_parts:
        return True, f"{len(new_parts)} new CC partition(s): {', '.join(new_parts)}"

    rmrb_size = os.path.getsize(rmrb_file) if rmrb_file and os.path.exists(rmrb_file) else 0
    if (
        os.path.basename(rmrb_file) != state.get("rmrb_file", "")
        or rmrb_size != state.get("rmrb_size", 0)
    ):
        return True, "RepBase combined file has changed"

    return False, "all partitions already merged with this RepBase version"


# ---------------------------------------------------------------------------
# Append logic (mirrors command_append from famdb.py)
# ---------------------------------------------------------------------------

def _parse_embl_version_from_header(rmrb_file):
    """Read the RELEASE line from the header before the first ID record."""
    with open(rmrb_file) as fh:
        for i, line in enumerate(fh):
            if line.startswith("ID"):
                break
            m = re.match(r"^(?:CC|##)\s+RELEASE\s+(\S+);", line)
            if m:
                return m.group(1)
            if i > 80:
                break
    return ""


def _parse_embl_copyright_header(rmrb_file):
    """Return the CC comment block that precedes the first ID record."""
    lines = []
    with open(rmrb_file) as fh:
        for line in fh:
            if line.startswith("ID"):
                break
            m = re.match(r"CC\s?(.*)", line)
            if m:
                text = m.group(1).rstrip("*").strip()
                lines.append(text)
    return "\n".join(lines)


def do_merge(db_dir, rmrb_file, dup_file, name=None, description=None):
    """
    Open the FamDB in read-write mode and append RepBase families from
    rmrb_file, skipping entries in dup_file and entries already present.

    Returns (added, total) counts.
    """
    # Load exclusion list (lowercase names to skip)
    exclusion_set = set()
    if dup_file and os.path.exists(dup_file):
        with open(dup_file) as fh:
            exclusion_set = {line.strip().lower() for line in fh if line.strip()}
        LOGGER.info(f"Loaded {len(exclusion_set):,} exclusion names from {dup_file}")
    else:
        LOGGER.warning("No exclusion file found; proceeding without it")

    LOGGER.info(f"Opening FamDB at {db_dir} for writing ...")
    db = FamDB(db_dir, "r+")

    lookup = db.get_all_taxa_names()
    LOGGER.info(f"  Taxonomy lookup: {len(lookup):,} entries")

    message = f"Adding Families From {os.path.basename(rmrb_file)}"
    changelog_rec = db.append_start_changelog(message)

    total_ctr = 0
    added_ctr = 0
    file_counts = {}
    new_val_taxa = set()
    dups = set()
    missing_parts = {}

    cc_components = db.components[COMPONENT_CC]

    for entry in FamDB.read_embl_families(rmrb_file, lookup):
        # Skip entries in the exclusion list or already present in any file
        if entry.accession.lower() in exclusion_set or not db.check_unique(entry):
            continue

        total_ctr += 1
        acc = entry.accession
        added = False

        # Route to every applicable CC partition that is locally installed
        add_leaves = {}
        add_taxa = set()
        for clade in entry.clades:
            part_dict = db.find_taxon(clade)
            cc_part = part_dict.get("cc") if part_dict else None
            if cc_part is not None and cc_part in cc_components:
                leaf = cc_components[cc_part]
                add_leaves[cc_part] = leaf
                if not db.get_families_for_taxon(clade):
                    add_taxa.add(clade)
            elif cc_part is not None:
                missing_parts[cc_part] = missing_parts.get(cc_part, 0) + 1

        if not add_leaves:
            LOGGER.debug(f"{acc}: no local CC partition found for its clade(s)")

        for part_num, leaf in add_leaves.items():
            try:
                leaf.add_family(entry)
                db.files[0]._add_family_taxon_links(acc, entry.clades)
                LOGGER.debug(f"Added {acc} to CC partition {part_num}")
                if not added:
                    added_ctr += 1
                    added = True
                file_counts[part_num] = file_counts.get(part_num, 0) + 1
            except Exception as exc:
                LOGGER.debug(f"Skipping duplicate {acc}: {exc}")
                dups.add(acc)

        if added:
            new_val_taxa.update(add_taxa)

    db.append_finish_changelog(message, changelog_rec)
    db.update_changelog(added_ctr, total_ctr, file_counts, rmrb_file)

    LOGGER.info(f"Added {added_ctr:,} of {total_ctr:,} families")
    if dups:
        LOGGER.debug(f"{len(dups)} duplicate accessions skipped")
    for part_num, cnt in missing_parts.items():
        LOGGER.info(
            f"CC partition {part_num} not installed locally; "
            f"{cnt} entries for that partition were not added"
        )

    # Update DB-level metadata
    db_info = db.get_metadata()
    if db_info:
        repbase_version = _parse_embl_version_from_header(rmrb_file)
        rb_copyright = _parse_embl_copyright_header(rmrb_file)

        new_name = name or db_info["name"]
        # Avoid appending duplicate copyright blocks
        new_copyright = db_info["copyright"]
        if rb_copyright and rb_copyright not in new_copyright:
            new_copyright = new_copyright + "\n\n" + rb_copyright

        new_desc = db_info["description"]
        if description:
            new_desc = new_desc + "\n" + description

        db.set_db_info(
            new_name,
            db_info["db_version"],
            db_info["date"],
            new_desc,
            new_copyright,
        )

    # Rebuild sparse taxonomy tree for nodes that gained their first family
    if new_val_taxa:
        LOGGER.info(f"Rebuilding pruned taxonomy tree ({len(new_val_taxa)} newly-valued nodes)")
        db.rebuild_pruned_tree(new_val_taxa)

    LOGGER.info("Finalizing files ...")
    db.finalize()
    db.close()

    return added_ctr, total_ctr


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--db-dir",
        help=(
            "Directory containing the FamDB files (*.0.h5 root + component files). "
            "Defaults to Libraries/famdb relative to the famdb.py installation directory."
        ),
    )
    parser.add_argument(
        "--meta",
        metavar="RMRBMeta.embl",
        help="Path to RMRBMeta.embl (metadata-only EMBL, included with RepeatMasker)",
    )
    parser.add_argument(
        "--seqs",
        metavar="RMRBSeqs.embl",
        help="Path to RMRBSeqs.embl (sequence EMBL, must be obtained from GIRI/RepBase)",
    )
    parser.add_argument(
        "--combined",
        metavar="RMRB.embl",
        help=(
            "Path to a pre-combined RMRB.embl (output of a previous --meta + --seqs run, "
            "or the file produced by RepeatMasker's addRepBase.pl). "
            "When provided, --meta and --seqs are ignored."
        ),
    )
    parser.add_argument(
        "--dup",
        metavar="RMRB_DUP.txt",
        help=(
            "Path to RMRB_DUP.txt exclusion list (families already curated into Dfam). "
            "Included with RepeatMasker in Libraries/RMRB_DUP.txt."
        ),
    )
    parser.add_argument(
        "--name",
        help="Override the database name (e.g. 'Dfam withRBRM')",
    )
    parser.add_argument(
        "--description",
        help="Additional text to append to the database description",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run the merge even if the state file says all partitions are up-to-date",
    )
    parser.add_argument(
        "-l", "--log-level",
        default="INFO",
        metavar="LEVEL",
        help="Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)",
    )
    return parser


def main():
    logging.basicConfig(format="%(levelname)s: %(message)s")
    parser = build_args()
    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper(), logging.INFO))

    if not args.db_dir:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        default = os.path.join(script_dir, "..", "Libraries", "famdb")
        if os.path.isdir(default):
            args.db_dir = default

    if not args.db_dir or not os.path.isdir(args.db_dir):
        LOGGER.error(
            "FamDB directory not found. Specify one with -i, or place the FamDB files "
            "in Libraries/famdb relative to the famdb.py installation directory."
        )
        sys.exit(1)

    db_dir = os.path.abspath(args.db_dir)

    # --- Sentinel check -------------------------------------------------
    sentinel = _sentinel_path(db_dir)
    if os.path.exists(sentinel):
        with open(sentinel) as fh:
            timestamp = fh.read().strip()
        LOGGER.error(
            f"A previous merge appears to have been interrupted on {timestamp}.\n"
            f"The files in {db_dir} may be in an inconsistent state.\n"
            "Remove all FamDB files, re-download the partitions, and re-run this script."
        )
        sys.exit(1)

    # Libraries/ is the canonical default location for all source files.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    libraries_dir = os.path.abspath(os.path.join(script_dir, "..", "Libraries"))
    LOGGER.info(f"Default Libraries directory: {libraries_dir}")

    def resolve_source(arg_path, filename, flag):
        """Return absolute path for a source file; log the origin or exit on error."""
        if arg_path:
            path = os.path.abspath(arg_path)
            if not os.path.isfile(path):
                LOGGER.error(f"File specified with {flag} not found: {path}")
                sys.exit(1)
            LOGGER.info(f"Sourcing {filename} from {flag} argument: {path}")
            return path
        default = os.path.join(libraries_dir, filename)
        if os.path.isfile(default):
            LOGGER.info(f"Sourcing {filename} from Libraries directory: {default}")
            return default
        return None

    # --- Locate RepBase files -------------------------------------------
    rmrb_file = None

    if args.combined:
        rmrb_file = resolve_source(args.combined, "RMRB.embl", "--combined")
        LOGGER.info(f"Using pre-combined RepBase file: {rmrb_file}")

    else:
        meta_file = resolve_source(args.meta, "RMRBMeta.embl", "--meta")
        seqs_file = resolve_source(args.seqs, "RMRBSeqs.embl", "--seqs")

        # Where to write / look for the combined file (alongside the source files)
        combined_path = os.path.join(
            os.path.dirname(meta_file) if meta_file else libraries_dir,
            "RMRB.embl",
        )

        if os.path.isfile(combined_path):
            LOGGER.info(f"Found existing combined file: {combined_path}")
            rmrb_file = combined_path
        elif meta_file and seqs_file:
            LOGGER.info(f"Combining {meta_file} + {seqs_file} -> {combined_path}")
            combine_rmrb_files(meta_file, seqs_file, combined_path)
            rmrb_file = combined_path
        elif meta_file and not seqs_file:
            LOGGER.error(
                "RMRBSeqs.embl not found in Libraries/ and --seqs was not supplied.\n"
                "Download RepBase RepeatMasker Edition from https://www.girinst.org/ "
                "and place it in Libraries/ or pass its path with --seqs.\n"
                "Alternatively, provide a pre-combined file with --combined."
            )
            sys.exit(1)
        else:
            LOGGER.error(
                "No RepBase EMBL files found in Libraries/ and no CLI paths were supplied.\n"
                "Provide --combined RMRB.embl, or both --meta RMRBMeta.embl and --seqs RMRBSeqs.embl."
            )
            sys.exit(1)

    # --- Locate exclusion list ------------------------------------------
    dup_file = resolve_source(args.dup, "RMRB_DUP.txt", "--dup")
    if not dup_file:
        LOGGER.warning(
            "RMRB_DUP.txt not found in Libraries/ and --dup was not supplied; "
            "proceeding without an exclusion list"
        )

    # --- State check ----------------------------------------------------
    state = load_state(db_dir)
    current_cc = get_cc_partition_files(db_dir)

    if not current_cc:
        LOGGER.warning("No curated-consensus (CC) component files found in %s", db_dir)
        LOGGER.warning("Download at least the CC partition files before running this script.")
        sys.exit(0)

    needs, reason = check_needs_merge(state, current_cc, rmrb_file)

    if not needs and not args.force:
        repbase_ver = state.get("repbase_version", "unknown")
        LOGGER.info(
            f"Nothing to do: {reason} (RepBase version: {repbase_ver}). "
            "Use --force to re-run anyway."
        )
        sys.exit(0)

    LOGGER.info(f"Merge needed: {reason}")

    # --- Write sentinel -------------------------------------------------
    with open(sentinel, "w") as fh:
        fh.write(str(datetime.datetime.now()) + "\n")

    # --- Perform merge --------------------------------------------------
    try:
        added, total = do_merge(
            db_dir,
            rmrb_file,
            dup_file,
            name=args.name,
            description=args.description,
        )
    except Exception:
        LOGGER.exception("Merge failed.  The sentinel file has been left in place.")
        sys.exit(1)

    # --- Update state ---------------------------------------------------
    repbase_version = _get_embl_version(rmrb_file)
    state["repbase_version"] = repbase_version
    state["rmrb_file"] = os.path.basename(rmrb_file)
    state["rmrb_size"] = os.path.getsize(rmrb_file)
    state["last_merge"] = str(datetime.datetime.now())
    # Record every CC file present at the time of the merge
    for fname, size in current_cc.items():
        state.setdefault("merged_cc", {})[fname] = {"size": size, "merged": True}

    save_state(db_dir, state)

    # --- Remove sentinel ------------------------------------------------
    os.remove(sentinel)

    LOGGER.info(
        f"Done. Added {added:,} RepBase families (RepBase version: {repbase_version or 'unknown'})."
    )


if __name__ == "__main__":
    main()
