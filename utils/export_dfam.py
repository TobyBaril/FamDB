#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Export the dfam database to FamDB format (v3).

    Usage: export_dfam.py [-h] [-l LOG_LEVEL]
               [--from-db mysql://...] [-r]
               --partition-dir DIR
               [--from-tax-dump ncbi_tax/]
               [--from-hmm file.hmm [--from-hmm file2.hmm ...]]
               [--db-version 3.2]
               [--db-date YYYY-MM-DD]
               --from-embl file.embl [--from-embl file2.embl ...]
               --repeat-peps file
               outfile

    In v3 format, families are split across component files by
    (curated/uncurated) x (consensus/hmm):

        outfile.0.h5                        -- root (taxonomy index only)
        outfile.curated.consensus.0.h5      -- all curated consensus seqs
        outfile.curated.hmm.N.h5            -- curated pHMMs, partition N
        outfile.uncurated.consensus.N.h5    -- uncurated consensus, partition N
        outfile.uncurated.hmm.N.h5          -- uncurated pHMMs, partition N

    The partition number N is determined by the partition_dfam.py script and
    the files it generates:

         - F_curated_consensus.json, F_curated_hmm.json
         - F_uncurated_consensus.json, F_uncurated_hmm.json

    Root File Creation (outfile.0.h5)

      Creates the index/navigation file that all clients open first. Contains:
      - Full taxonomy tree -- all NCBI taxonomy nodes with parent/child links,
        scientific names, and aliases for Dfam-related taxa (those appearing in
        any family's clades or their ancestors).
      - RepeatPeps -- the protein sequences from RepeatPeps.lib, used for
        protein-based repeat masking
      - File map -- metadata describing all the component files that make up
        this FamDB release (partition assignments, partition roots, file names)
      - Database metadata -- version, date, copyright, description

      No family data lives here -- it's pure navigation and reference data.

    The export runs in three phases:

      Phase 1 -- Pre-computation (main process):
        Queries the family_clade join table twice (once for DF/curated, once
        for DR/uncurated), yielding (family_id, taxon_id) pairs. Each pair is
        routed to the appropriate component and partition using the F.json maps.
        Result:
           a dict of {component -> {partition_num -> [family_id, ...]}}.

      Phase 2 -- Parallel pickle workers:
        One worker process per partition chunk. Each worker:
        1. Connects to the DB independently
        2. Fetches full family records by primary key (consensus sequence or HMM
           profile + all annotations/metadata)
        3. For HMM workers: stores the gzip-compressed HMM blob as-is (defers
           decompression to Phase 3 to keep pickle files small)
        4. Writes a stream of serialized Family objects to a temp .pkl file

        Data produced: temporary .pkl files on disk, one per chunk. NOTE: The
        HMM is stored in the database's compressed state.

      Phase 3 -- Sequential HDF5 assembly:
        One process per partition, reads its pickle chunk(s) and writes the
        final .h5 leaf file. Each leaf file contains:
          - File metadata -- partition key, component type, DB version/date,
            file map
          - Families dataset -- each family stored as an HDF5 group with
            datasets/attributes for:
            - Accession, name, description, classification
            - Consensus sequence (consensus files) or decompressed HMM profile
              (HMM files)
            - Clade associations (list of taxon IDs)
            - Citations, aliases, search stages, gathering thresholds, etc.
          - Pruned taxonomy tree -- a sparse version of the taxonomy containing
            only nodes relevant to this partition's families (for fast lineage
            lookups)
          - File changelog -- timestamped history of operations applied to this
            file

          Also ingests any families from --from-embl and --from-hmm flat files
          that belong to this partition's node set.

    Run partition_dfam.py once before this script -- it writes all four
    partition files to a single directory (default: "partitions"):
        ./partition_dfam.py -o partitions

    This generates:
        partitions/F_curated_hmm.json
        partitions/F_curated_consensus.json
        partitions/F_uncurated_hmm.json
        partitions/F_uncurated_consensus.json

    Pass that same directory here via --partition-dir.

    Data source options:

    --partition-dir     : Directory written by partition_dfam.py (required)
    -p, --partition          : Specify partition numbers to export (default: all)
    --from-db                : Connection string to MySQL database
    -t, --test-set           : Only include 5 families per partition
    --from-tax-dump          : NCBI taxonomy dump directory
    --from-embl              : EMBL-format file(s) to import
    --from-hmm               : HMM-format file(s) to import
    --repeat-peps            : Path to RepeatPeps.lib
    -c --dfam-config         : Path to dfam config file
    --db-version             : Database version string
    --db-date                : Database date (YYYY-MM-DD)
    --workers                : Max parallel pickle workers (default: cpu_count).
                               Lower values reduce memory at the cost of runtime.
    --chunk-size             : Max families per pickle worker (default: 500000).
                               Large partitions are split into chunks of this size.

SEE ALSO:
    famdb.py
    Dfam: http://www.dfam.org
"""

import argparse
import datetime
import itertools
import json
import logging
import pickle
import re
import tempfile
import time
import sys
import os
import shutil
import h5py
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

# Import SQLAlchemy
from sqlalchemy import create_engine, select, func
from sqlalchemy.orm import sessionmaker

# Import FamDB
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from famdb_classes import FamDB, FamDBLeaf, FamDBRoot
LOGGER = logging.getLogger(__name__)

from famdb_globals import (
    COPYRIGHT_TEXT,
    DESCRIPTION,
    META_META,
    META_DB_VERSION,
    META_DB_DATE,
    META_FILE_MAP,
    META_FILE_INFO,
    META_DB_NAME,
    META_DB_COPYRIGHT,
    META_DB_DESCRIPTION,
    COMPONENT_CC,
    COMPONENT_CH,
    COMPONENT_UC,
    COMPONENT_UH,
    COMPONENT_META,
    COMPONENT_TYPES,
    GROUP_FAMILIES,
)
from famdb_helper_methods import families_iterator, accession_bin


# Dfam admin dependencies -- set PYTHONPATH to include Dfam Schemata/ORMs/python and Lib
try:
    import dfamorm as dfam
    import DfamConfig as dc
except ImportError as e:
    sys.exit(
        f"Missing Dfam admin dependency: {e}\n"
        "This utility requires internal Dfam libraries.\n"
        "Add the Dfam Schemata/ORMs/python and Lib directories to PYTHONPATH."
    )

from famdb_data_loaders import (
    load_taxonomy_from_db,
    load_taxonomy_from_dump,
    iterate_db_families,
    iterate_db_families_by_ids,
    read_hmm_families,
)

# Maps component type to (curated_str, model_str) used in filenames
COMPONENT_FILE_PARTS = {
    COMPONENT_CC: ("curated",   "consensus"),
    COMPONENT_CH: ("curated",   "hmm"),
    COMPONENT_UC: ("uncurated", "consensus"),
    COMPONENT_UH: ("uncurated", "hmm"),
}


def _taxon_scientific_name(tax_db, tax_id):
    """Return the scientific name for a tax_id from tax_db, or str(tax_id)."""
    node = tax_db.get(tax_id)
    if node:
        for name_class, name_val in node.names:
            if name_class == "scientific name":
                return name_val
    return str(tax_id)


def build_component_file_map(out_str, cc_F, ch_F, uc_F, uh_F, tax_db):
    """
    Build the META_FILE_MAP dict for all v3 component files.

    Keys follow the pattern "<component_type>.<partition_num>" (e.g. "cc.0",
    "ch.1") plus "0" for the root file.
    """
    file_map = {}

    # Root entry
    root_tax_id = 1
    file_map["0"] = {
        "T_root": root_tax_id,
        "filename": f"{os.path.basename(out_str)}.0.h5",
        "F_roots": [root_tax_id],
        "T_root_name": _taxon_scientific_name(tax_db, root_tax_id),
        "F_roots_names": [],
    }

    component_F_map = [
        (COMPONENT_CC, cc_F),
        (COMPONENT_CH, ch_F),
        (COMPONENT_UC, uc_F),
        (COMPONENT_UH, uh_F),
    ]
    for comp_type, F_dict in component_F_map:
        if not F_dict:
            continue
        curated_str, model_str = COMPONENT_FILE_PARTS[comp_type]
        for n, partition in F_dict.items():
            key = f"{comp_type}.{n}"
            t_root = partition["T_root"]
            f_roots = partition.get("F_roots", [])
            f_roots_names = [
                _taxon_scientific_name(tax_db, r)
                for r in f_roots
                if r != t_root
            ]
            file_map[key] = {
                "T_root": t_root,
                "filename": f"{os.path.basename(out_str)}.{curated_str}.{model_str}.{n}.h5",
                "F_roots": f_roots,
                "T_root_name": _taxon_scientific_name(tax_db, t_root),
                "F_roots_names": f_roots_names,
            }
    return file_map


def build_partition_cache(cc_nodes, ch_nodes, uc_nodes, uh_nodes):
    """
    Build the PartitionCache dict from per-component partition node dicts.

    Each nodes_dict maps {partition_num: [tax_id, ...]} for one component.
    Returns {tax_id_str: {"cc": N|None, "ch": N|None, "uc": N|None, "uh": N|None}}.
    """
    def invert(nodes_dict):
        result = {}
        for part_num, node_list in nodes_dict.items():
            for node in node_list:
                result[str(node)] = int(part_num)
        return result

    cc_map = invert(cc_nodes)
    ch_map = invert(ch_nodes)
    uc_map = invert(uc_nodes)
    uh_map = invert(uh_nodes)

    all_nodes = set(cc_map) | set(ch_map) | set(uc_map) | set(uh_map)
    cache = {}
    for node in all_nodes:
        entry = {
            COMPONENT_CC: cc_map.get(node),
            COMPONENT_CH: ch_map.get(node),
            COMPONENT_UC: uc_map.get(node),
            COMPONENT_UH: uh_map.get(node),
        }
        if any(v is not None for v in entry.values()):
            cache[node] = entry
    return cache


def collect_family_taxon_map(filename):
    """
    Scan a component HDF5 file and return a {tax_id_str: [accession, ...]} dict.
    Used in the post-processing step to build the root's Lookup/ByTaxon.
    """
    mapping = {}
    with h5py.File(filename, "r") as f:
        if GROUP_FAMILIES not in f:
            return mapping
        for acc in families_iterator(f[GROUP_FAMILIES], GROUP_FAMILIES):
            path = accession_bin(acc)
            entry = f.get(f"{path}/{acc}")
            if entry is not None:
                clades = entry.attrs.get("clades")
                if clades is not None:
                    for clade_id in clades:
                        mapping.setdefault(str(int(clade_id)), []).append(acc)
    return mapping


def precompute_family_partitions(session, cc_nodes, ch_nodes, uc_nodes, uh_nodes,
                                 test_limit=None):
    """
    Pre-compute which families belong to which component+partition.

    Performs one bulk JOIN query for curated (DF%) families and one for
    uncurated (DR%) families, then routes each (family_id, tax_id) pair to
    the appropriate component+partition using the inverted F.json maps.

    Returns {component_type: {partition_num: sorted_list_of_family_ids}}.
    A family may appear in multiple partitions if its clades span partition
    boundaries (matching the existing per-partition query semantics).
    """
    def invert(nodes_dict):
        """Map {tax_id: partition_num} from {partition_num: [tax_id]}."""
        result = {}
        if not nodes_dict:
            return result
        for part_num, node_list in nodes_dict.items():
            for tax_id in node_list:
                result[tax_id] = int(part_num)
        return result

    cc_inv = invert(cc_nodes)
    ch_inv = invert(ch_nodes)
    uc_inv = invert(uc_nodes)
    uh_inv = invert(uh_nodes)

    # Accumulate family IDs per (component, partition) using sets for dedup
    assignments = {comp: defaultdict(set) for comp in COMPONENT_TYPES}

    def _query_and_assign(accession_prefix, inv_maps):
        """Run one JOIN query and distribute results into assignments."""
        labelled = [(comp, inv) for comp, inv in inv_maps if inv]
        if not labelled:
            return
        LOGGER.info(
            "Pre-computing %s family->partition assignments", accession_prefix
        )
        start = time.perf_counter()
        count = 0
        query = (
            select(
                dfam.Family.id,
                dfam.t_family_clade.c.dfam_taxdb_tax_id,
            )
            .join(
                dfam.t_family_clade,
                dfam.Family.id == dfam.t_family_clade.c.family_id,
            )
            .where(dfam.Family.accession.like(f"{accession_prefix}%"))
            .where(dfam.Family.disabled != 1)
        )
        if test_limit:
            query = query.limit(test_limit * 500)  # rough upper bound
        for row in session.execute(query).yield_per(10000):
            count += 1
            fid, tax_id = row.id, row.dfam_taxdb_tax_id
            for comp, inv in labelled:
                part = inv.get(tax_id)
                if part is not None:
                    assignments[comp][part].add(fid)
        delta = time.perf_counter() - start
        LOGGER.info(
            "Assigned %d (family, clade) rows for %s in %.1fs",
            count, accession_prefix, delta,
        )

    _query_and_assign("DF", [(COMPONENT_CC, cc_inv), (COMPONENT_CH, ch_inv)])
    _query_and_assign("DR", [(COMPONENT_UC, uc_inv), (COMPONENT_UH, uh_inv)])

    # Convert sets to sorted lists; optionally trim for test mode
    result = {}
    for comp, part_map in assignments.items():
        result[comp] = {}
        for part, id_set in part_map.items():
            ids = sorted(id_set)
            if test_limit:
                ids = ids[:test_limit]
            result[comp][part] = ids

    return result


# ---------------------------------------------------------------------------
# Phase 2: chunk pickle workers (must be module-level to be picklable)
# ---------------------------------------------------------------------------

def run_chunk_worker(db_conn_str, family_ids, is_hmm, pickle_path, label, batch_size):
    """
    Fetch families by primary key and write pickled Family objects to pickle_path.

    This is a top-level function so ProcessPoolExecutor can pickle it.
    Each worker connects to the DB independently, queries by ID (no complex
    JOIN), and streams results directly to a pickle file.
    """
    engine = create_engine(db_conn_str)
    session = sessionmaker(engine)()

    count = 0
    start = time.perf_counter()
    report_every = max(1000, len(family_ids) // 100)
    t_pickle_write = 0.0

    try:
        # For HMM components, keep the model as compressed bytes in the pickle
        # so workers skip decompression CPU cost.  The compressed blob is
        # stored directly in HDF5; get_family() decompresses on read.
        defer_decompress = is_hmm

        LOGGER.debug("%s: Fetching %d families", label, len(family_ids))
        with open(pickle_path, "wb") as pf:
            for item in iterate_db_families_by_ids(
                session, family_ids, is_hmm=is_hmm, batch_size=batch_size,
                defer_model_decompress=defer_decompress,
            ):
                _pw = time.perf_counter()
                pickle.dump(item, pf)
                t_pickle_write += time.perf_counter() - _pw
                count += 1
                if count % report_every == 0:
                    elapsed = time.perf_counter() - start
                    avg_ms = elapsed / count * 1000
                    rem = int(avg_ms / 1000 * (len(family_ids) - count))
                    LOGGER.debug(
                        "%s: %5d / %d : %dms avg : %s remaining",
                        label, count, len(family_ids), int(avg_ms),
                        str(datetime.timedelta(seconds=rem)),
                    )
        total_elapsed = time.perf_counter() - start
        t_db = total_elapsed - t_pickle_write
        LOGGER.debug(
            "%s: Pickled %d families in %s  [db/build=%.1fs  pickle_write=%.1fs]",
            label, count,
            str(datetime.timedelta(seconds=int(total_elapsed))),
            t_db, t_pickle_write,
        )
        return {"families": count, "db_build_time": t_db, "pickle_write_time": t_pickle_write, "elapsed": total_elapsed}
    finally:
        session.close()
        engine.dispose()


# ---------------------------------------------------------------------------
# Phase 3: assemble HDF5 from pickle chunks
# ---------------------------------------------------------------------------

# Fields whose content is the bulk sequence/model data -- excluded from the
# metadata size tally so we measure only the annotation overhead.
_META_SIZE_SKIP = frozenset({"consensus", "model"})


def _metadata_field_sizes(family):
    """Return {field_name: estimated_bytes} for set metadata fields, excluding
    consensus and model.  Uses vars() to read directly from __dict__ so it
    works even when family.model is None (deferred decompress path).

    Size estimates:
      str   -> UTF-8 encoded byte length
      list  -> 8 bytes x element count (covers clades: list of ints)
      int / float / bool -> 8 bytes (fixed-width scalar)
    """
    sizes = {}
    for k, v in vars(family).items():
        if k in _META_SIZE_SKIP or v is None:
            continue
        if isinstance(v, str):
            sizes[k] = len(v.encode("utf-8"))
        elif isinstance(v, list):
            sizes[k] = 8 * len(v)
        elif isinstance(v, (int, float, bool)):
            sizes[k] = 8
    return sizes


def build_component_from_pickles(
    args, component_type, partition_num, chunk_files,
    out_str, nodes, metadata, tax_lookup
):
    """
    Stream through sorted pickle chunk files and write a single HDF5 leaf file.

    Also appends any --from-embl and --from-hmm families that belong to this
    partition's node set, matching the behaviour of the old export_families().
    """
    curated_str, model_str = COMPONENT_FILE_PARTS[component_type]
    filename = f"{out_str}.{curated_str}.{model_str}.{partition_num}.h5"
    label = os.path.basename(filename)
    LOGGER.debug("Building %s from %d chunk(s)", label, len(chunk_files))

    part_key = f"{component_type}.{partition_num}"
    start = time.perf_counter()
    count = 0

    with FamDBLeaf(filename, "w", component_type=component_type) as outfile:
        outfile.set_metadata(
            part_key,
            metadata[META_FILE_INFO],
            metadata[META_DB_NAME],
            metadata[META_DB_VERSION],
            metadata[META_DB_DATE],
            metadata[META_DB_COPYRIGHT],
        )

        message = "Preparing To Add Family Data"
        log_timestamp = outfile.update_changelog(message)

        # Stream pickle chunks in order.
        # HMM items are (Family, raw_gzip_bytes) tuples; the compressed blob is
        # passed directly to add_family() which stores it as-is in HDF5.
        # No decompression happens here -- get_family() decompresses on read.
        t_pickle_read = 0.0
        t_hdf5_write = 0.0
        meta_field_totals = defaultdict(int)   # {field_name: total_bytes}
        hmm_compressed_bytes = 0
        for chunk_path in sorted(chunk_files):
            with open(chunk_path, "rb") as pf:
                while True:
                    try:
                        _pr = time.perf_counter()
                        item = pickle.load(pf)
                        t_pickle_read += time.perf_counter() - _pr

                        if isinstance(item, tuple):
                            family, blob = item
                            if blob:
                                family.model = blob  # keep as compressed bytes
                                hmm_compressed_bytes += len(blob)
                        else:
                            family = item

                        for field, sz in _metadata_field_sizes(family).items():
                            meta_field_totals[field] += sz

                        _hw = time.perf_counter()
                        outfile.add_family(family)
                        t_hdf5_write += time.perf_counter() - _hw
                        count += 1
                    except EOFError:
                        break

        # EMBL file imports (typically /dev/null in production)
        for embl_file in args.from_embl:
            for family in FamDB.read_embl_families(embl_file, tax_lookup, nodes):
                outfile.add_family(family)
                count += 1

        # HMM file imports
        for hmm_file in args.from_hmm:
            for family in read_hmm_families(hmm_file, tax_lookup, nodes):
                outfile.add_family(family)
                count += 1

        delta = str(datetime.timedelta(seconds=int(time.perf_counter() - start)))
        outfile._verify_change(log_timestamp, message)
        outfile.update_changelog(f"Added {count} Families In {delta}", verified=True)
        outfile.finalize()

    total_elapsed = time.perf_counter() - start
    delta = str(datetime.timedelta(seconds=int(total_elapsed)))
    hdf5_other = total_elapsed - t_pickle_read - t_hdf5_write
    LOGGER.debug(
        "Written: %s (%d families in %s)  "
        "[pickle_read=%.1fs  hdf5_write=%.1fs  other=%.1fs"
        "  hdf5=%.2fms/family]",
        label, count, delta,
        t_pickle_read, t_hdf5_write, hdf5_other,
        t_hdf5_write / count * 1000 if count else 0,
    )
    if count and meta_field_totals:
        total_meta = sum(meta_field_totals.values())
        avg_bytes = total_meta / count
        top = sorted(meta_field_totals.items(), key=lambda x: -x[1])
        field_summary = "  ".join(f"{k}={v/count:.0f}B" for k, v in top)
        LOGGER.debug(
            "Metadata size %s: %.0f B/family avg (%.1f MB total)  %s",
            label, avg_bytes, total_meta / 1e6, field_summary,
        )
    if hmm_compressed_bytes:
        LOGGER.debug(
            "HMM compressed model storage %s: %.1f MB total  (%.0f B/family)",
            label,
            hmm_compressed_bytes / 1e6,
            hmm_compressed_bytes / count,
        )
    return {
        "families": count,
        "pickle_read_time": t_pickle_read,
        "hdf5_write_time": t_hdf5_write,
        "elapsed": total_elapsed,
        "meta_field_totals": dict(meta_field_totals),
        "hmm_compressed_bytes": hmm_compressed_bytes,
    }


def file_check(infile):
    try:
        with open(infile):
            pass
        return True
    except IOError:
        return False


def run_export_root(args, out_str, tax_db, metadata, conf):
    """Create the root file with full taxonomy and RepeatPeps (no families)."""
    LOGGER.info("Exporting root file")
    root_path = f"{out_str}.0.h5"
    with FamDBRoot(root_path, "w") as root:
        root.write_full_taxonomy(tax_db)
        root.write_repeatpeps(args.repeat_peps)
        root.set_metadata(
            "0",
            metadata[META_FILE_INFO],
            metadata[META_DB_NAME],
            metadata[META_DB_VERSION],
            metadata[META_DB_DATE],
            metadata[META_DB_COPYRIGHT],
        )
        root.finalize()
    LOGGER.info(f"Root file written: {root_path}")


def main():
    """Parses command-line arguments and runs the import."""

    logging.basicConfig()
    _main_start = time.perf_counter()

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--log-level", default="INFO")
    parser.add_argument("--partition-dir", required=True,
                        help="Directory written by partition_dfam.py "
                             "(must contain F_curated_hmm.json etc.)")
    parser.add_argument("-p", "--partition", nargs="+", default=[])
    parser.add_argument("-t", "--test-set", action="store_true", default=False)
    parser.add_argument("--from-tax-dump")
    parser.add_argument("--from-embl", action="append", default=[], required=True)
    parser.add_argument("--from-hmm", action="append", default=[])
    parser.add_argument("-c", "--dfam-config", dest="dfam_config")
    parser.add_argument("--repeat-peps", dest="repeat_peps", required=True)
    parser.add_argument("--db-version")
    parser.add_argument("--db-date")
    parser.add_argument("--count-taxa-in")
    parser.add_argument("--min-init", action="store_true")
    parser.add_argument(
        "--components",
        nargs="+",
        choices=["cc", "ch", "uc", "uh"],
        default=None,
        metavar="COMP",
        help=(
            "Limit which component files are built.  Accepts one or more of "
            "'cc' (curated consensus), 'ch' (curated HMM), 'uc' (uncurated consensus), "
            "'uh' (uncurated HMM).  The root index file is always written.  "
            "Default: all four components."
        ),
    )
    parser.add_argument(
        "--root-only",
        action="store_true",
        default=False,
        help=(
            "Write only the root index file (taxonomy, RepeatPeps, empty "
            "Lookup/ByTaxon and PartitionCache) -- no component family files.  "
            "Equivalent to --min-init without the sequencing-artifacts rename."
        ),
    )
    parser.add_argument("--workers", type=int, default=None,
                        help="Max parallel pickle workers (default: cpu_count). "
                             "Reduce to limit memory usage.")
    parser.add_argument("--chunk-size", type=int, default=500_000,
                        help="Max families per pickle worker chunk (default: 500000). "
                             "Large partitions are split into chunks of this size.")
    parser.add_argument("--batch-size", type=int, default=500,
                        help="Families per SQL batch within each worker (default: 500).")
    parser.add_argument("-o", "--outfile", required=True,
                        help="Output file base name (e.g. 'dfam' produces dfam.0.h5, dfam.curated.hmm.0.h5, ...)")

    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    if not os.path.isdir(args.partition_dir):
        LOGGER.error(
            f"Partition directory '{args.partition_dir}' does not exist. "
            "Run partition_dfam.py first."
        )
        exit(1)

    needed_files = (
        [f for f in [args.dfam_config, args.repeat_peps] if f]
        + args.from_embl
        + args.from_hmm
    )
    for file in needed_files:
        if not file_check(file):
            LOGGER.error(f"Could not read {file}, export canceled")
            exit(1)

    conf = dc.DfamConfig(args.dfam_config)

    df_engine = create_engine(conf.getDBConnStrWPassFallback("Dfam"))
    df_sfactory = sessionmaker(df_engine)
    session = df_sfactory()

    db_url = df_engine.url
    LOGGER.info(
        f"Connected to schema '{db_url.database}' on {db_url.host}"
        + (f":{db_url.port}" if db_url.port else "")
        + f" as '{db_url.username}'"
    )

    db_conn_str = conf.getDBConnStrWPassFallback("Dfam")

    db_version = None
    db_date = None

    version_info = session.execute(select(dfam.DbVersion)).scalar_one()
    db_version = version_info.dfam_version
    db_date = version_info.dfam_release_date.strftime("%Y-%m-%d")

    if args.db_version:
        db_version = args.db_version
    if args.db_date:
        db_date = args.db_date

    if not db_version:
        raise Exception("Could not determine database version.")
    if not db_date:
        db_date = datetime.date.today().strftime("%Y-%m-%d")

    year_match = re.match(r"(\d{4})-", db_date)
    if year_match:
        db_year = year_match.group(1)
    else:
        raise Exception("Date should be in YYYY-MM-DD format, got: " + db_date)

    copyright_text = COPYRIGHT_TEXT % (db_year, db_version, db_date)
    out_str = args.outfile

    # --- Load partition F files from --partition-dir ---
    def load_F(label):
        path = os.path.join(args.partition_dir, f"F_{label}.json")
        if not os.path.exists(path):
            LOGGER.warning(
                f"Partition file not found: {path} -- "
                f"skipping {label} component."
            )
            return None, None, None
        with open(path, "r") as fh:
            data = json.load(fh)
        F = data["F"]
        F_meta = data[META_META]
        if F_meta[META_DB_VERSION] != db_version or F_meta[META_DB_DATE] != db_date:
            LOGGER.error(
                f"Partition file {path} does not match current database version "
                f"(file: {F_meta[META_DB_VERSION]} / {F_meta[META_DB_DATE]}, "
                f"db: {db_version} / {db_date}). "
                f"Re-run partition_dfam.py before export."
            )
            sys.exit(1)
        return F, F_meta, {int(k): v["nodes"] for k, v in F.items()}

    # Determine which components are active.
    # --root-only or --components not given with no value -> no components.
    # --components cc ch ... -> explicit subset.
    # Neither flag -> all four.
    if args.root_only:
        active_components = set()
    elif args.components is not None:
        active_components = set(args.components)
    else:
        active_components = set(COMPONENT_TYPES)   # all four

    if active_components != set(COMPONENT_TYPES):
        LOGGER.info(
            "Selective export -- active components: %s",
            ", ".join(sorted(active_components)) if active_components else "(none -- root only)",
        )

    _no = (None, None, None)  # sentinel for skipped component load
    cc_F, cc_meta, cc_nodes = load_F("curated_consensus")   if COMPONENT_CC in active_components else _no
    ch_F, ch_meta, ch_nodes = load_F("curated_hmm")         if COMPONENT_CH in active_components else _no
    uc_F, uc_meta, uc_nodes = load_F("uncurated_consensus") if COMPONENT_UC in active_components else _no
    uh_F, uh_meta, uh_nodes = load_F("uncurated_hmm")       if COMPONENT_UH in active_components else _no

    if active_components and not any([cc_F, ch_F, uc_F, uh_F]):
        LOGGER.error(
            f"No partition files found in '{args.partition_dir}'. "
            "Run partition_dfam.py first."
        )
        exit(1)

    # Collect all relevant taxonomy nodes from all four partition files
    relevant_nodes = []
    for nodes_dict in [cc_nodes, ch_nodes, uc_nodes, uh_nodes]:
        if nodes_dict:
            for node_list in nodes_dict.values():
                relevant_nodes.extend(node_list)
    relevant_nodes = list(set(relevant_nodes))

    if not args.from_tax_dump:
        tax_db, tax_lookup = load_taxonomy_from_db(session, relevant_nodes)
    else:
        tax_db, tax_lookup = load_taxonomy_from_dump(args.from_tax_dump, relevant_nodes)

    # Build unified file map for all component files
    file_map = build_component_file_map(
        out_str,
        cc_F or {},
        ch_F or {},
        uc_F or {},
        uh_F or {},
        tax_db,
    )
    F_meta = cc_meta or ch_meta or uc_meta or uh_meta
    file_info = {META_META: F_meta, META_FILE_MAP: file_map}

    meta_database = "Dfam"
    description = DESCRIPTION
    if args.min_init:
        meta_database = "Sequencing_artifacts_only"
        db_version = "1.0"
        description = "A minimal library composed of sequencing artifacts only."

    metadata = {
        META_FILE_INFO: file_info,
        META_DB_NAME: meta_database,
        META_DB_VERSION: db_version,
        META_DB_DATE: db_date,
        META_DB_DESCRIPTION: description,
        META_DB_COPYRIGHT: copyright_text,
    }

    # --- Export root file (no families) ---
    run_export_root(args, out_str, tax_db, metadata, conf)

    if args.min_init:
        elapsed = time.perf_counter() - _main_start
        LOGGER.info("--min-init: skipping family export -- done in %s",
                    str(datetime.timedelta(seconds=int(elapsed))))
        return

    if not active_components:
        elapsed = time.perf_counter() - _main_start
        LOGGER.info("Root-only export complete -- done in %s",
                    str(datetime.timedelta(seconds=int(elapsed))))
        return

    # --- Phase 1: Pre-compute family->partition assignments ---
    test_limit = 5 if args.test_set else None
    family_partitions = precompute_family_partitions(
        session,
        cc_nodes or {},
        ch_nodes or {},
        uc_nodes or {},
        uh_nodes or {},
        test_limit=test_limit,
    )

    # Filter to requested partition numbers if -p was given
    if args.partition:
        requested = set(int(p) for p in args.partition)
        for comp in list(family_partitions):
            family_partitions[comp] = {
                p: ids for p, ids in family_partitions[comp].items()
                if p in requested
            }

    # Log partition sizes
    for comp in COMPONENT_TYPES:
        for part, ids in sorted(family_partitions[comp].items()):
            curated_str, model_str = COMPONENT_FILE_PARTS[comp]
            LOGGER.info(
                "Partition %s.%s.%d: %d families",
                curated_str, model_str, part, len(ids),
            )

    # Phase 1 is complete -- release the DB connection now so it doesn't time
    # out during the long Phase 2/3 workers and cause a spurious error on exit.
    session.close()
    df_engine.dispose()

    # --- Phase 2: Build chunk jobs and run pickle workers ---
    tmp_dir = tempfile.mkdtemp(prefix=f"{os.path.basename(out_str)}_chunks_",
                               dir=os.path.dirname(out_str) or ".")
    LOGGER.info("Temp chunk directory: %s", tmp_dir)

    # chunk_map: {(comp, part_num): [pickle_path, ...]}
    chunk_map = defaultdict(list)

    # Build chunk jobs: consensus first so they claim worker slots immediately
    def _make_chunk_jobs(comp_list):
        jobs = []
        _, is_hmm_flag = None, None
        for comp in comp_list:
            is_curated, is_hmm_flag = COMPONENT_META[comp]
            curated_str, model_str = COMPONENT_FILE_PARTS[comp]
            for part_num, family_ids in sorted(family_partitions[comp].items()):
                if not family_ids:
                    continue
                chunks = [
                    family_ids[i:i + args.chunk_size]
                    for i in range(0, len(family_ids), args.chunk_size)
                ]
                for chunk_idx, chunk_ids in enumerate(chunks):
                    pkl = os.path.join(
                        tmp_dir,
                        f"{curated_str}.{model_str}.{part_num}.{chunk_idx:04d}.pkl",
                    )
                    chunk_map[(comp, part_num)].append(pkl)
                    label = (
                        f"{os.path.basename(out_str)}"
                        f".{curated_str}.{model_str}.{part_num}"
                        + (f" chunk{chunk_idx}" if len(chunks) > 1 else "")
                    )
                    jobs.append((comp, db_conn_str, chunk_ids, is_hmm_flag, pkl, label))
        return jobs

    consensus_jobs = _make_chunk_jobs([c for c in [COMPONENT_CC, COMPONENT_UC] if c in active_components])
    hmm_jobs       = _make_chunk_jobs([c for c in [COMPONENT_CH, COMPONENT_UH] if c in active_components])
    all_chunk_jobs = consensus_jobs + hmm_jobs

    all_complete = True
    executor = ProcessPoolExecutor(max_workers=args.workers)
    futures = {}
    phase2_start = time.perf_counter()
    families_total_p2 = sum(len(fids) for _, _, fids, _, _, _ in all_chunk_jobs)
    families_done_p2 = 0
    jobs_total_p2 = len(all_chunk_jobs)
    jobs_done_p2 = 0
    component_chunk_stats = defaultdict(list)
    try:
        for comp, db_conn, fids, is_hmm, pkl, label in all_chunk_jobs:
            future = executor.submit(
                run_chunk_worker, db_conn, fids, is_hmm, pkl, label, args.batch_size
            )
            futures[future] = (label, len(fids), comp)

        for future in as_completed(futures):
            label, n_fam, comp = futures[future]
            try:
                stats = future.result()
                families_done_p2 += n_fam
                jobs_done_p2 += 1
                if stats:
                    component_chunk_stats[comp].append(stats)
                p2_elapsed = time.perf_counter() - phase2_start
                avg_ms = p2_elapsed / families_done_p2 * 1000 if families_done_p2 else 0
                rem_fam = families_total_p2 - families_done_p2
                est_rem = datetime.timedelta(seconds=int(avg_ms / 1000 * rem_fam))
                LOGGER.info(
                    "Phase 2 (%d/%d chunks): %s/%s families done -- %.1fms avg -- %s remaining",
                    jobs_done_p2, jobs_total_p2,
                    f"{families_done_p2:,}", f"{families_total_p2:,}",
                    avg_ms, est_rem,
                )
            except Exception:
                all_complete = False
                LOGGER.error("Error in chunk worker %s:\n%s", label, traceback.format_exc())

    except KeyboardInterrupt:
        LOGGER.warning("Interrupted -- cancelling pending jobs and terminating workers...")
        executor.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        executor.shutdown(wait=False)

    for comp in COMPONENT_TYPES:
        stats_list = component_chunk_stats.get(comp, [])
        if not stats_list:
            continue
        total_fam = sum(s["families"] for s in stats_list)
        if not total_fam:
            continue
        total_db = sum(s["db_build_time"] for s in stats_list)
        total_pw  = sum(s["pickle_write_time"] for s in stats_list)
        curated_str, model_str = COMPONENT_FILE_PARTS[comp]
        LOGGER.info(
            "Phase 2 summary -- %s.%s: %d chunks, %s families -- %.1fms avg/family "
            "(db_build=%.1fms  pickle_write=%.1fms)",
            curated_str, model_str, len(stats_list), f"{total_fam:,}",
            (total_db + total_pw) / total_fam * 1000,
            total_db / total_fam * 1000,
            total_pw / total_fam * 1000,
        )

    if not all_complete:
        LOGGER.error("Chunk worker errors encountered. HDF5 assembly skipped.")
        LOGGER.info("Temp files left for inspection: %s", tmp_dir)
        return

    # --- Phase 3: Assemble HDF5 files from pickle chunks (parallel) ---
    LOGGER.info("Assembling HDF5 component files from pickle chunks")

    assembly_jobs = []
    nodes_by_comp = {
        COMPONENT_CC: cc_nodes, COMPONENT_CH: ch_nodes,
        COMPONENT_UC: uc_nodes, COMPONENT_UH: uh_nodes,
    }
    for comp in COMPONENT_TYPES:
        if comp not in active_components:
            continue
        for part_num, family_ids in sorted(family_partitions[comp].items()):
            if not family_ids:
                continue
            chunks = chunk_map.get((comp, part_num), [])
            if not chunks:
                LOGGER.warning("No chunk files found for %s.%d -- skipping", comp, part_num)
                continue
            nodes_dict = nodes_by_comp[comp]
            nodes_for_part = nodes_dict.get(part_num, []) if nodes_dict else []
            assembly_jobs.append((comp, part_num, chunks, nodes_for_part))

    all_assembly_complete = True
    assembly_executor = ProcessPoolExecutor(max_workers=args.workers)
    assembly_futures = {}
    phase3_start = time.perf_counter()
    families_total_p3 = sum(len(family_partitions[comp][part_num]) for comp, part_num, _, _ in assembly_jobs)
    families_done_p3 = 0
    jobs_total_p3 = len(assembly_jobs)
    jobs_done_p3 = 0
    component_assembly_stats = defaultdict(list)
    try:
        for comp, part_num, chunks, nodes_for_part in assembly_jobs:
            future = assembly_executor.submit(
                build_component_from_pickles,
                args, comp, part_num, chunks,
                out_str, nodes_for_part, metadata, tax_lookup,
            )
            assembly_futures[future] = (comp, part_num, len(family_partitions[comp][part_num]))

        for future in as_completed(assembly_futures):
            comp, part_num, n_fam = assembly_futures[future]
            label = f"{comp}.{part_num}"
            try:
                stats = future.result()
                families_done_p3 += n_fam
                jobs_done_p3 += 1
                if stats:
                    component_assembly_stats[comp].append(stats)
                p3_elapsed = time.perf_counter() - phase3_start
                avg_ms = p3_elapsed / families_done_p3 * 1000 if families_done_p3 else 0
                rem_fam = families_total_p3 - families_done_p3
                est_rem = datetime.timedelta(seconds=int(avg_ms / 1000 * rem_fam))
                LOGGER.info(
                    "Phase 3 (%d/%d partitions): %s/%s families done -- %.1fms avg -- %s remaining",
                    jobs_done_p3, jobs_total_p3,
                    f"{families_done_p3:,}", f"{families_total_p3:,}",
                    avg_ms, est_rem,
                )
            except Exception:
                all_assembly_complete = False
                LOGGER.error("Error assembling %s:\n%s", label, traceback.format_exc())

    except KeyboardInterrupt:
        LOGGER.warning("Interrupted -- cancelling assembly jobs...")
        assembly_executor.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        assembly_executor.shutdown(wait=False)

    for comp in COMPONENT_TYPES:
        stats_list = component_assembly_stats.get(comp, [])
        if not stats_list:
            continue
        total_fam = sum(s["families"] for s in stats_list)
        if not total_fam:
            continue
        total_hdf5 = sum(s["hdf5_write_time"] for s in stats_list)
        agg_meta = defaultdict(int)
        for s in stats_list:
            for field, sz in s["meta_field_totals"].items():
                agg_meta[field] += sz
        hmm_bytes = sum(s["hmm_compressed_bytes"] for s in stats_list)
        curated_str, model_str = COMPONENT_FILE_PARTS[comp]
        meta_total = sum(agg_meta.values())
        top_fields = sorted(agg_meta.items(), key=lambda x: -x[1])
        field_summary = "  ".join(f"{k}={v/total_fam:.0f}B" for k, v in top_fields)
        LOGGER.info(
            "Phase 3 summary -- %s.%s: %d partitions, %s families -- "
            "hdf5=%.1fms/family avg -- metadata=%.0f B/family avg  %s",
            curated_str, model_str, len(stats_list), f"{total_fam:,}",
            total_hdf5 / total_fam * 1000,
            meta_total / total_fam if total_fam else 0,
            field_summary,
        )
        if hmm_bytes:
            LOGGER.info(
                "Phase 3 summary -- %s.%s: HMM compressed storage %.1f MB total (%.0f B/family avg)",
                curated_str, model_str, hmm_bytes / 1e6, hmm_bytes / total_fam,
            )

    if not all_assembly_complete:
        LOGGER.error("Assembly errors encountered. Post-processing skipped.")
        LOGGER.info("Pickle chunks left for inspection: %s", tmp_dir)
        return

    # Clean up temp directory
    try:
        shutil.rmtree(tmp_dir)
        LOGGER.info("Removed temp chunk directory: %s", tmp_dir)
    except Exception as exc:
        LOGGER.warning("Could not remove temp dir %s: %s", tmp_dir, exc)

    # --- Post-processing: write Lookup/ByTaxon and PartitionCache to root ---
    LOGGER.info("Post-processing: building Lookup/ByTaxon and PartitionCache")

    # Build the list of component files that were actually created
    component_files = []
    for comp in COMPONENT_TYPES:
        curated_str, model_str = COMPONENT_FILE_PARTS[comp]
        for part_num in sorted(family_partitions[comp].keys()):
            component_files.append((comp, part_num))

    family_taxon_map = {}
    for comp, part_num in component_files:
        curated_str, model_str = COMPONENT_FILE_PARTS[comp]
        filename = f"{out_str}.{curated_str}.{model_str}.{part_num}.h5"
        if os.path.exists(filename):
            LOGGER.info(f"Scanning {filename} for family-taxon map")
            for tax_id_str, accessions in collect_family_taxon_map(filename).items():
                for acc in accessions:
                    family_taxon_map.setdefault(tax_id_str, [])
                    if acc not in family_taxon_map[tax_id_str]:
                        family_taxon_map[tax_id_str].append(acc)

    # Build PartitionCache from each component's node dicts
    partition_cache = build_partition_cache(
        {int(k): v for k, v in cc_nodes.items()} if cc_nodes else {},
        {int(k): v for k, v in ch_nodes.items()} if ch_nodes else {},
        {int(k): v for k, v in uc_nodes.items()} if uc_nodes else {},
        {int(k): v for k, v in uh_nodes.items()} if uh_nodes else {},
    )

    with FamDBRoot(f"{out_str}.0.h5", "r+") as root:
        root.write_lookup_bytaxon(family_taxon_map)
        root.write_partition_cache(partition_cache)

    # --- Move all files to output directory ---
    famdb_dir = f"./{out_str}"
    LOGGER.info(f"Moving files to {famdb_dir}")
    os.makedirs(famdb_dir, exist_ok=True)
    for file in os.listdir("."):
        if file.startswith(os.path.basename(out_str)) and file.endswith(".h5"):
            shutil.move(f"./{file}", f"{famdb_dir}/{file}")

    # --- Build pruned tree ---
    LOGGER.info(f"Building pruned tree in {famdb_dir}")
    new_famdb = FamDB(famdb_dir, "r+")
    new_famdb.build_pruned_tree()

    new_famdb.print_info(history=True)

    elapsed = time.perf_counter() - _main_start
    LOGGER.info("Export complete.")
    LOGGER.info("Run time: %s", datetime.timedelta(seconds=elapsed))


if __name__ == "__main__":
    main()
