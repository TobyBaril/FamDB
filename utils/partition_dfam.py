#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Usage: ./partition_dfam.py [--help] [--log-level] [--dfam-config] [--version]
                               [--output-dir DIR] [--chunk-size #] [--rep-base]
                               [--root-def] [--save-csv] [--save-tree]

    This script partitions the Dfam database into taxonomy groups of approximately
    equal file size.  All output is written to a single directory (--output-dir,
    default: "partitions").  Pass that same directory to export_dfam.py via
    --partition-dir and it will find everything it needs automatically.

    OUTPUT FILES
    ============

    Always written (required by export_dfam.py):

        <output-dir>/F_curated_hmm.json
        <output-dir>/F_curated_consensus.json
        <output-dir>/F_uncurated_hmm.json
        <output-dir>/F_uncurated_consensus.json

            Partition assignment for each (curation x model) combination.
            Each file contains a metadata block and an "F" dict mapping
            partition number -> {T_root, bytes, nodes, F_roots}.
            "nodes" is the list of NCBI tax_ids whose families belong in
            that partition file; "F_roots" is the topmost node of the
            partition's subtree (entry point from the taxonomy tree).

    Written by --save-csv:

        <output-dir>/T_orig.csv

            Per-node table of unweighted filesizes (all zeros before weights
            are applied).  Useful for verifying the topology was built correctly.

        <output-dir>/T_{label}_partitioned.csv

            Per-node table (node, chunk, filesize) after partition assignment
            for each combination.  Useful for diagnosing imbalanced partitions
            or verifying that a taxon landed in the expected chunk.

    Written by --save-tree (requires --save-csv for the CSV input to R):

        <output-dir>/T_{label}.newick

            Newick-format representation of the full taxonomy tree with
            partition chunk ids as node labels.  Can be loaded into any
            tree viewer (e.g. FigTree, iTOL) for visual inspection.

        <output-dir>/T_{label}.png

            PNG visualization of the newick tree coloured by chunk,
            produced by tree.R via Rscript.  Only written when Rscript
            is on PATH and --save-csv is also given (the R script reads
            the partitioned CSV for colouring).

    Intermediate caches (always written, safe to delete between runs):

        <output-dir>/topo.pkl           -- topology (shared across all four combinations)
        <output-dir>/nodes_{label}.pkl  -- per-taxon byte estimates per combination
        <output-dir>/RMRB_sizes.json    -- RepBase taxon sizes (only with --rep_base)

    ALGORITHM
    =========

    Node weights differ per (curation_status, model_type) combination:

        curated     : only DF* families
        uncurated   : only DR* families
        consensus   : weight = SUM(family.length)  [sequence bytes]
        hmm         : weight = SUM(length + length*177 + 1160 + len(description))
                       [estimated HMM blob bytes]

    A single bottom-up O(N) pass accumulates subtree weight from leaves to root.
    When a node's accumulated weight first reaches --chunk_size, a partition
    boundary is cut there and its weight is not propagated further up the tree.
    A final top-down pass labels every node with its assigned partition number.

    Args:
        --help, -h        : Show this help message and exit
        --log-level, -l   : Control the logger level of the script
        --dfam-config, -c : Dfam Config file
        --version, -v     : Get Dfam Version
        --output-dir, -o  : Directory to write all output files (default: "partitions")
        --chunk-size, -S  : Target partition size in bytes (default 70,000,000,000 = 70 GB)
        --rep-base, -r    : Path to RMRB.embl; reserves space for RepBase data in partitions
        --root-def, -d    : File listing tax_ids (one per line) to pre-assign to partition 0
        --save-csv        : Write per-node CSV tables (T_orig.csv, T_*_partitioned.csv)
        --save-tree       : Write newick tree and PNG visualization (requires --save-csv for PNG colouring)

SEE ALSO: related_script.py
          Dfam: http://www.dfam.org

AUTHOR(S):
    Anthony Gray agray@systemsbiology.org

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
# Module imports
import argparse
import copy
import logging
import pickle
import os
import sys
import json
import time
import uuid
import subprocess
import shutil

# SQL Alchemy
from sqlalchemy import create_engine, text, select
from sqlalchemy.orm import sessionmaker

# Dfam admin dependencies -- set PYTHONPATH to include Dfam Schemata/ORMs/python and Lib
try:
    import dfamorm as dfam
    import DfamConfig as dc
    import DfamVersion as dfVersion
except ImportError as e:
    sys.exit(
        f"Missing Dfam admin dependency: {e}\n"
        "This utility requires internal Dfam libraries.\n"
        "Add the Dfam Schemata/ORMs/python and Lib directories to PYTHONPATH."
    )

# Import FamDB
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from famdb_globals import META_DB_DATE, META_DB_VERSION, META_META, META_UUID

LOGGER = logging.getLogger(__name__)

# All (curation_status, model_type) combinations to partition
COMBINATIONS = [
    ("curated",   "hmm"),
    ("curated",   "consensus"),
    ("uncurated", "hmm"),
    ("uncurated", "consensus"),
]


def _fmt_elapsed(seconds):
    h = int(seconds) // 3600
    m = (int(seconds) % 3600) // 60
    s = seconds % 60
    if h:
        return f"{h}:{m:02d}:{s:05.2f}"
    if m:
        return f"{m}:{s:05.2f}"
    return f"{s:.2f}s"


def _usage():
    """Print out docstring as program usage"""
    help(os.path.splitext(os.path.basename(__file__))[0])
    sys.exit(0)


def parse_RMRB(args, session, RB_file):
    data = []
    with open(args.rep_base, "r") as input:
        lines = input.readlines()
        fam = {"species": None, "seq_size": 0}
        for line in lines:
            if line.startswith("CC        Species:"):
                fam["species"] = [
                    spec.strip() for spec in line.split(":")[1].split(",")
                ]
            elif line.startswith("SQ   Sequence"):
                fam["seq_size"] = int(line.split(" ")[4]) * 8

            if fam["species"] and fam["seq_size"]:
                data.append(fam)
                fam = {"species": None, "seq_size": 0}

    looked_up = {}
    with session.bind.begin() as conn:
        for fam in data:
            fam["tax_id"] = []
            for species in fam["species"]:
                if species in looked_up:
                    tax_id = looked_up[species]
                else:
                    query = f"SELECT tax_id FROM `ncbi_taxdb_names` WHERE sanitized_name='{species}'"
                    res = tuple(conn.execute(text(query)))
                    tax_id = int(res[0][0]) if res else None
                    looked_up[species] = tax_id
                    if not tax_id:
                        LOGGER.warning(f"Could not resolve RepBase species '{species}' to an NCBI taxon; skipping")
                fam["tax_id"] += [tax_id]
    with open(RB_file, "w+") as output:
        output.write(json.dumps(data))


def build_T_topology(args, session, topo_file, RB_file):
    """
    Build (and cache) the taxonomy tree topology from the Dfam + NCBI databases.

    The returned T dict has filesize=0 and tot_weight=0 for all nodes;
    only the tree structure (parent, children) is established here.
    This is called once and shared across all four (curation x model) combinations.
    """
    if os.path.exists(topo_file):
        LOGGER.info("Found Stashed Topology")
        with open(topo_file, "rb") as phandle:
            return pickle.load(phandle)

    LOGGER.info("Building Tree Topology")

    node_query = (
        "SELECT dfam_taxdb.tax_id, parent_id FROM `ncbi_taxdb_nodes` "
        "JOIN dfam_taxdb ON dfam_taxdb.tax_id = ncbi_taxdb_nodes.tax_id"
    )

    if args.rep_base:
        with open(RB_file, "rb") as spec_file:
            RB_json = json.load(spec_file)
            RB_taxa = set(
                str(tid)
                for i in RB_json
                for tid in i["tax_id"]
                if tid is not None
            )
        node_query += (
            f" UNION SELECT tax_id, parent_id from ncbi_taxdb_nodes"
            f" WHERE tax_id IN ({','.join(RB_taxa)})"
        )

    with session.bind.begin() as conn:
        tax_ids, parent_ids = zip(*conn.execute(text(node_query)))
        tax_ids = [int(x) for x in tax_ids]
        parent_ids = [int(x) for x in parent_ids]

        while True:
            missing_parents = [p for p in parent_ids if p not in tax_ids]
            if not missing_parents:
                break
            update_query = (
                f"SELECT tax_id, parent_id FROM `ncbi_taxdb_nodes`"
                f" WHERE tax_id IN ({','.join(str(n) for n in missing_parents)})"
            )
            new_taxs, new_parents = zip(*conn.execute(text(update_query)))
            tax_ids.extend([int(x) for x in new_taxs])
            parent_ids.extend([int(x) for x in new_parents])

    T = {
        z[0]: {
            "parent": z[1],
            "children": [],
            "filesize": 0,
            "tot_weight": 0,
            "chunk": None,
        }
        for z in zip(tax_ids, parent_ids)
    }

    for n in T:
        parent = T[n]["parent"]
        if parent:
            T[parent]["children"].append(n)

    T[1]["children"].remove(1)
    T[1]["parent"] = None

    LOGGER.info("Stashing Topology")
    with open(topo_file, "wb") as phandle:
        pickle.dump(T, phandle, protocol=4)
    return T


def query_filesizes(session, model_type, curation_status, node_file):
    """
    Query per-taxon byte-size estimates for a given (model_type, curation_status)
    combination.  Results are cached to node_file.

    curation_status : "curated" (DF* only), "uncurated" (DR* only)
    model_type      : "consensus" (sequence length) or "hmm" (blob size estimate)
    """
    if os.path.exists(node_file):
        LOGGER.debug(f"Found stashed node sizes ({curation_status}/{model_type})")
        with open(node_file, "rb") as phandle:
            return pickle.load(phandle)  # (filesizes, famcounts) tuple

    if model_type == "consensus":
        size_expr = "SUM(COALESCE(family.length, 0))"
    else:
        # COALESCE(family.length, 0) is required: uncurated families often have
        # NULL length, and NULL arithmetic propagates to NULL, causing SUM() to
        # silently ignore those families and massively underestimate partition 0.
        size_expr = (
            "SUM((COALESCE(family.length, 0) * 178 + 1160"
            " + OCTET_LENGTH(COALESCE(family.description, ''))))"
        )

    if curation_status == "curated":
        curation_filter = "AND family.accession LIKE 'DF%'"
    else:
        curation_filter = "AND family.accession LIKE 'DR%'"

    node_query = (
        f"SELECT family_clade.dfam_taxdb_tax_id,"
        f" {size_expr} AS byte_est,"
        f" COUNT(DISTINCT family.id) AS fam_count"
        f" FROM family_clade JOIN family ON family_clade.family_id = family.id"
        f" WHERE family.disabled != 1 {curation_filter}"
        f" GROUP BY family_clade.dfam_taxdb_tax_id"
    )

    LOGGER.debug(f"Querying node sizes ({curation_status}/{model_type})")
    with session.bind.begin() as conn:
        filesizes = {}
        famcounts = {}
        for row in conn.execute(text(node_query)):
            filesizes[int(row[0])] = int(row[1]) if row[1] else 0
            famcounts[int(row[0])] = int(row[2]) if row[2] else 0

        # Check once per curation_status (hmm is queried first in COMBINATIONS).
        # Families with NULL or zero length are silently zero-weighted; warn so
        # the operator knows the size estimates may not reflect export reality.
        if model_type == "hmm":
            null_query = (
                f"SELECT COUNT(*) FROM family_clade"
                f" JOIN family ON family_clade.family_id = family.id"
                f" WHERE family.disabled != 1 {curation_filter}"
                f" AND (family.length IS NULL OR family.length = 0)"
            )
            (null_count,) = conn.execute(text(null_query)).one()
            LOGGER.debug(
                f"  null/zero-length check ({curation_status}): {null_count:,} families"
            )
            if null_count:
                LOGGER.warning(
                    f"{null_count:,} {curation_status} families have NULL or zero length "
                    f"-- their size contribution is estimated as HMM overhead only ({1160} bytes each). "
                    f"Partition size estimates may be inaccurate for nodes that hold these families."
                )

    with open(node_file, "wb") as phandle:
        pickle.dump((filesizes, famcounts), phandle, protocol=4)
    return filesizes, famcounts


def _post_order(T, root):
    """
    Iterative post-order traversal of T starting at root.
    Yields each node after all of its descendants have been yielded.
    Avoids Python recursion limits for large taxonomy trees.
    """
    stack = [(root, False)]
    while stack:
        n, visited = stack.pop()
        if visited:
            yield n
        else:
            stack.append((n, True))
            for child in reversed(T[n]["children"]):
                stack.append((child, False))


def apply_weights(T_topo, filesizes, famcounts, RB_json=None):
    """
    Return a fresh weighted copy of T_topo with filesize, famcount, and tot_weight populated.

    T_topo   : base topology dict (not modified)
    filesizes: {tax_id: bytes} from query_filesizes
    famcounts: {tax_id: int} family count per taxon, from query_filesizes
    RB_json  : parsed RMRB_sizes.json list, or None
    """
    T = copy.deepcopy(T_topo)

    for n in T:
        T[n]["famcount"] = 0

    for taxon, size in filesizes.items():
        if taxon in T:
            T[taxon]["filesize"] = size

    for taxon, count in famcounts.items():
        if taxon in T:
            T[taxon]["famcount"] = count

    if RB_json:
        for fam in RB_json:
            for tax_id in fam["tax_id"]:
                if tax_id is None:
                    continue
                if tax_id in T:
                    T[tax_id]["filesize"] += fam["seq_size"]
                else:
                    T[tax_id]["filesize"] = fam["seq_size"]

    # Iterative post-order accumulation (avoids recursion limit on deep trees)
    for n in _post_order(T, 1):
        T[n]["tot_weight"] = T[n]["filesize"] + sum(
            T[c]["tot_weight"] for c in T[n]["children"]
        )

    return T


def run_chunk_assignment(T, S, args):
    """
    Bottom-up O(N) chunk assignment.

    A single post-order pass accumulates subtree weight upward; a cut is made
    whenever a node's accumulated weight meets the size threshold S (and the
    node is not the tree root or pre-assigned to chunk 0 via --root_def).
    A final top-down DFS then labels every node with its chunk id.

    Returns F: {chunk_id: {T_root, bytes, nodes, F_roots}}.
    T is modified in-place (T[n]["chunk"] is set for every n).
    """

    # --- Handle --root_def pre-assignments (chunk 0) ---
    pre_assigned = set()
    root_offset = 0

    if args.root_def:
        with open(args.root_def, "r") as f:
            root_ids = [int(line.strip()) for line in f if line.strip()]

        # Mark entire subtrees as pre-assigned to chunk 0 (iterative)
        for rid in root_ids:
            stack = [rid]
            while stack:
                n = stack.pop()
                pre_assigned.add(n)
                stack.extend(T[n]["children"])
            root_offset += T[rid]["tot_weight"]

        LOGGER.info(
            f"Chunk 0 pre-assigned weight: {root_offset / 1_073_741_824:,.2f} Gb"
        )

    # --- Single post-order pass: accumulate weight, cut when threshold met ---
    # accumulated[n] = total uncut weight rooted at n (own filesize + uncut children)
    accumulated = {n: (0 if n in pre_assigned else T[n]["filesize"]) for n in T}
    chunk_cut = {}   # {node: chunk_id} -- nodes where a partition boundary is cut
    chunk_ctr = 1

    for n in _post_order(T, 1):
        if n not in pre_assigned and n != 1 and accumulated[n] >= S:
            # Cut a new chunk rooted here
            chunk_cut[n] = chunk_ctr
            LOGGER.debug(
                f"  cut chunk {chunk_ctr}: node {n}, "
                f"size {accumulated[n] / 1_073_741_824:,.2f} GB"
            )
            chunk_ctr += 1
            accumulated[n] = 0          # consumed by this cut; do not propagate up

        parent = T[n]["parent"]
        if parent is not None:
            accumulated[parent] += accumulated[n]

    LOGGER.debug(
        f"  chunk 0 remainder: {(accumulated[1] + root_offset) / 1_073_741_824:,.2f} GB"
    )

    # --- Top-down DFS: assign chunk labels to every node ---
    stack = [(1, 0)]   # (node, inherited_chunk_id)
    while stack:
        n, inherited = stack.pop()
        if n in pre_assigned:
            my_chunk = 0
        elif n in chunk_cut:
            my_chunk = chunk_cut[n]
        else:
            my_chunk = inherited
        T[n]["chunk"] = my_chunk
        for child in T[n]["children"]:
            stack.append((child, my_chunk))

    # --- Build F dict ---
    F = {
        cid: {"T_root": 1 if cid == 0 else None, "bytes": 0, "nodes": [], "F_roots": []}
        for cid in range(chunk_ctr)
    }
    for n, cid in chunk_cut.items():
        F[cid]["T_root"] = n

    for n in T:
        cid = T[n]["chunk"]
        F[cid]["nodes"].append(n)
        F[cid]["bytes"] += T[n]["filesize"]

    # --- Build F_roots ---
    # Chunk 0 always includes the tree root (node 1)
    F[0]["F_roots"].append(1)

    # Non-zero chunks: the F_root is always the cut node -- the topmost node of
    # that chunk's contiguous subtree.  (In bottom-up partitioning, a cut
    # node's parent may itself be inside another cut chunk, so we can't rely
    # on "parent is in chunk 0" as the old top-down algorithm did.)
    for n, cid in chunk_cut.items():
        F[cid]["F_roots"].append(n)

    return F


def save_combination_outputs(
    T, F, db_version, db_date, curation_status, model_type,
    output_dir, session, args
):
    """
    Save the F JSON, partitioned CSV, and newick/PNG for one combination.
    Also logs a summary of chunk sizes and roots.
    """
    label = f"{curation_status}_{model_type}"
    LOGGER.info(f"Partitions ({label}):")

    for n in sorted(F.keys()):
        size_gb = sum(T[i]["filesize"] for i in F[n]["nodes"]) / 1_073_741_824
        tax_names = []
        for i in F[n]["F_roots"]:
            tax_rec = session.execute(
                select(dfam.DfamTaxdb).where(dfam.DfamTaxdb.tax_id == i)
            ).scalar_one_or_none()
            tax_name = ""
            if tax_rec:
                tax_name = tax_rec.scientific_name
                if tax_rec.common_name is not None:
                    tax_name = tax_name + " (" + tax_rec.common_name + ")"
                tax_name = tax_name + " [" + str(i) + "]"
            else:
                tax_rec = session.execute(
                    select(dfam.NcbiTaxdbNames)
                    .where(dfam.NcbiTaxdbNames.name_class == "scientific name")
                    .where(dfam.NcbiTaxdbNames.tax_id == i)
                ).scalar_one_or_none()
                if tax_rec:
                    tax_name = tax_rec.name_txt
                tax_rec = (
                    session.execute(
                        select(dfam.NcbiTaxdbNames)
                        .where(dfam.NcbiTaxdbNames.name_class == "common name")
                        .where(dfam.NcbiTaxdbNames.tax_id == i)
                    )
                    .scalars()
                    .first()
                )
                if tax_rec:
                    tax_name = tax_name + " (" + tax_rec.name_txt + ")"
                tax_name = tax_name + " [" + str(i) + "]"
            if tax_name == "":
                tax_name = str(i)
            tax_names.append(tax_name)
        total_fams = sum(T[i]["famcount"] for i in F[n]["nodes"])
        root_str = ", ".join(tax_names) if tax_names else "(none)"
        LOGGER.info(f"  [{n:3d}]  {size_gb:8.2f} GB  {total_fams:>10,} families  {root_str}")

    # Save F JSON
    F_file = {
        META_META: {
            META_UUID: str(uuid.uuid4()),
            META_DB_VERSION: db_version,
            META_DB_DATE: db_date,
        },
        "F": F,
    }
    f_path = os.path.join(output_dir, f"F_{label}.json")
    with open(f_path, "w") as outfile:
        json.dump(F_file, outfile)
    LOGGER.info(f"  Saved {f_path}")

    if args.save_csv:
        csv_path = os.path.join(output_dir, f"T_{label}_partitioned.csv")
        with open(csv_path, "w") as outfile:
            outfile.write(
                "node, chunk, filesize\n"
                + "\n".join([f"{n},{T[n]['chunk']},{T[n]['filesize']}" for n in T])
            )
        LOGGER.info(f"  Saved {csv_path}")

    if args.save_tree:
        # Build newick (iterative post-order to avoid recursion limit on large trees)
        def build_newick(root):
            buf = {}
            for n in _post_order(T, root):
                children = T[n]["children"]
                if children:
                    buf[n] = "(" + ",".join(buf.pop(c) for c in children) + f"){n}"
                else:
                    buf[n] = str(n)
            return buf[root]

        newick_path = os.path.join(output_dir, f"T_{label}.newick")
        with open(newick_path, "w") as outfile:
            outfile.write(build_newick(1) + ";\n")
        LOGGER.info(f"  Saved {newick_path}")

        if args.save_csv:
            rscript_path = shutil.which("Rscript")
            if rscript_path is None:
                LOGGER.warning(
                    "Rscript not found on PATH - skipping PNG visualization. "
                    "Install R + ggtree if you want PNG output."
                )
            else:
                try:
                    png_path = os.path.join(output_dir, f"T_{label}.png")
                    subprocess.run(
                        [rscript_path, "tree.R", newick_path, csv_path, png_path],
                        check=True,
                    )
                    LOGGER.info(f"  Saved {png_path}")
                except subprocess.CalledProcessError as e:
                    LOGGER.error(f"Rscript execution failed: {e}")
        else:
            LOGGER.info("  Skipping PNG (--save-csv required for R colouring data)")


#
# main subroutine ( protected from import execution )
#
def main(*args):
    """Parse arguments and run"""

    logging.basicConfig(stream=sys.stdout, format="%(levelname)s: %(message)s")

    class _CustomUsageAction(argparse.Action):
        def __init__(
            self, option_strings, dest, default=False, required=False, help=None
        ):
            super(_CustomUsageAction, self).__init__(
                option_strings=option_strings,
                dest=dest,
                nargs=0,
                const=True,
                default=default,
                required=required,
                help=help,
            )

        def __call__(self, parser, args, values, option_string=None):
            _usage()

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action=_CustomUsageAction)
    parser.add_argument("-l", "--log-level", default="INFO")
    parser.add_argument("-c", "--dfam-config")
    parser.add_argument("-v", "--version", dest="get_version", action="store_true")
    parser.add_argument("-o", "--output-dir", default="partitions",
                        help="Directory to write all output files (default: partitions)")
    parser.add_argument("-S", "--chunk-size", default=70000000000)
    parser.add_argument("-r", "--rep-base")
    parser.add_argument("-d", "--root-def")
    parser.add_argument("--save-csv", action="store_true",
                        help="Write per-node CSV tables to output dir")
    parser.add_argument("--save-tree", action="store_true",
                        help="Write newick tree and PNG visualization to output dir")
    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    os.makedirs(args.output_dir, exist_ok=True)

    conf = dc.DfamConfig(args.dfam_config)
    df_ver = dfVersion.DfamVersion()
    version = df_ver.version_string

    if args.get_version:
        LOGGER.info(version)
        exit(0)

    dfamdb = create_engine(conf.getDBConnStrWPassFallback("Dfam"))
    dfamdb_sfactory = sessionmaker(dfamdb)
    session = dfamdb_sfactory()

    db_url = dfamdb.url
    LOGGER.info(
        f"Connected to schema '{db_url.database}' on {db_url.host}"
        + (f":{db_url.port}" if db_url.port else "")
        + f" as '{db_url.username}'"
    )

    version_info = session.execute(select(dfam.DbVersion)).scalar_one()
    db_version = version_info.dfam_version
    db_date = version_info.dfam_release_date.strftime("%Y-%m-%d")

    # ~ PARSE RMRB.embl ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RB_file = os.path.join(args.output_dir, "RMRB_sizes.json")
    if args.rep_base:
        if os.path.exists(RB_file):
            LOGGER.info("Found RepBase File")
        else:
            LOGGER.info("Generating RMRB.sizes")
            parse_RMRB(args, session, RB_file)

    RB_json = None
    if args.rep_base and os.path.exists(RB_file):
        with open(RB_file, "rb") as f:
            RB_json = json.load(f)

    # ~ BUILD TOPOLOGY (once, shared across all combinations) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    topo_file = os.path.join(args.output_dir, "topo.pkl")
    T_topo = build_T_topology(args, session, topo_file, RB_file)

    if args.save_csv:
        orig_csv_path = os.path.join(args.output_dir, "T_orig.csv")
        with open(orig_csv_path, "w") as outfile:
            outfile.write(
                "node, weight\n"
                + "\n".join([f"{n},{T_topo[n]['tot_weight']}" for n in T_topo])
            )
        LOGGER.info(f"Saved {orig_csv_path}")

    S = int(args.chunk_size)
    total_start = time.perf_counter()

    # ~ RUN ALL FOUR COMBINATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for curation_status, model_type in COMBINATIONS:
        label = f"{curation_status}_{model_type}"
        combo_start = time.perf_counter()
        LOGGER.info(f"\n=== Partitioning: {label} ===")

        node_file = os.path.join(args.output_dir, f"nodes_{label}.pkl")
        filesizes, famcounts = query_filesizes(session, model_type, curation_status, node_file)

        T = apply_weights(T_topo, filesizes, famcounts, RB_json)

        total_gb = T[1]["tot_weight"] / 1_073_741_824
        threshold_gb = S / 1_073_741_824
        LOGGER.info(f"Total size: {total_gb:,.2f} GB  (threshold: {threshold_gb:.2f} GB)")

        F = run_chunk_assignment(T, S, args)

        n_partitions = len(F)
        if n_partitions == 1:
            LOGGER.info("  ->1 partition (total size below threshold)")
        else:
            LOGGER.info(f"  ->{n_partitions} partitions")

        save_combination_outputs(
            T, F, db_version, db_date, curation_status, model_type,
            args.output_dir, session, args
        )

        elapsed = time.perf_counter() - combo_start
        LOGGER.info(f"  Completed in {_fmt_elapsed(elapsed)}")

    total_elapsed = time.perf_counter() - total_start
    LOGGER.info(f"\nAll combinations complete.  Total runtime: {_fmt_elapsed(total_elapsed)}")


#
# Wrap script functionality in main() to avoid automatic execution
# when imported ( e.g. when help is called on file )
#
if __name__ == "__main__":
    main(*sys.argv)
