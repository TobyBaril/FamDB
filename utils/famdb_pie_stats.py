#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    Usage: ./famdb_pie_stats.py [--help] [--log-level LEVEL] [--dfam-config CONFIG]
                                <command> [command-options]

    Generate the taxonomy distribution JSON used by the Dfam website pie chart,
    or discover a new set of high-level groups from the live database.

    Commands:
        generate    Compute per-group family counts and write the pie-stats JSON.
        suggest     Discover candidate group lists from the live database.

    GENERATE COMMAND
    ================
    Counts distinct TE families (DF* and DR* separately) with at least one
    family_clade entry within each named group's subtree.  Families with no clade
    in any named group are reported under "Other".  A family may appear in more
    than one group if it spans multiple major clades.  Groups with count == 0 are
    omitted.

    When --count-by species is given, counts distinct species-rank taxa
    (ncbi_taxdb_nodes.rank = 'species') that have at least one family_clade entry
    within the group's subtree, rather than counting the families themselves.

    When --groups-file is given, the curated and uncurated taxon lists are read
    from that JSON file instead of the hard-coded CURATED_GROUPS / UNCURATED_GROUPS
    constants.  The file format is identical to this script's output (or the output
    of the suggest command), so a previous run can be fed back in to regenerate
    counts after a database update.  Existing count values and any "Other" entry
    (taxon 0) in the file are ignored and recalculated automatically.

    Output format:
        {
          "curated":   [{"group": "...", "taxon": <id>, "count": <n>}, ...],
          "uncurated": [{"group": "...", "taxon": <id>, "count": <n>}, ...],
          "meta": {"db_version": "...", "db_date": "...", "count_by": "families"|"species"}
        }

    SUGGEST COMMAND
    ===============
    Discovers candidate group lists independently for curated (DF*) and uncurated
    (DR*) families using a top-down greedy cut on the taxonomy tree.

    A node becomes a group boundary when its subtree contains >= --min-fraction
    of that dataset's families AND no individual child subtree is also that large
    (i.e., it is the deepest node that is still large enough as a coherent unit).

    Counts used for discovery are approximate (per-node family-clade sums; a
    family with N clade entries is counted N times across ancestors).  The final
    JSON output always uses accurate DISTINCT counts.

    Prints a human-readable table for both datasets, then writes the JSON using
    the discovered groups (suitable for feeding back into "generate --groups-file").

    Args:
        --help, -h           : Show this help message and exit
        --log-level, -l      : Logger level (default: INFO)
        --dfam-config, -c    : Dfam Config file

    generate args:
        --groups-file, -g    : JSON file with curated/uncurated taxon lists
        --output, -o         : Output JSON file (default: pie_stats.json)
        --count-by, -b       : Count "families" (default) or "species"

    suggest args:
        --output, -o         : Output JSON file (default: pie_stats.json)
        --min-fraction, -f   : Minimum fraction of families for a group (default: 0.03)
        --max-groups, -m     : Maximum number of groups to suggest (default: 20)
        --count-by, -b       : Count "families" (default) or "species" in the final output

SEE ALSO: partition_dfam.py
          Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley rhubley@systemsbiology.org

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
import argparse
import json
import logging
import os
import sys

from sqlalchemy import create_engine, text, select
from sqlalchemy.orm import sessionmaker

LOGGER = logging.getLogger(__name__)

# Groups for the curated (DF*) pie chart.
# Run with --suggest to discover whether a new set fits the current database better.
# Groups with count == 0 are omitted.  "Other" (taxon 0) is appended automatically.
CURATED_GROUPS = [
    ("Mammals",              40674),
    ("Ray-Finned Fish",      7898),
    ("Butterflies & Moths",  104431),
    ("Beetles",              7041),
    ("Flies",                7147),
    ("Roundworms",           6231),
    ("Birds",                8782),
]

# Groups for the uncurated (DR*) pie chart.
UNCURATED_GROUPS = [
    ("Mammals",              40674),
    ("Birds",                8782),
    ("Lepidosaurs",          8504),
    ("Turtles",              8459),
    ("Ray-Finned Fish",      7898),
    ("Amphibians",           8292),
    ("Plants",               33090),
    ("Butterflies & Moths",  104431),
    ("Ants, Bees & Wasps",   7399),
    ("Flies",                7147),
]

# Maximum IN-clause length before chunking (keeps queries within MySQL packet limits)
_CHUNK = 8000


# ---------------------------------------------------------------------------
# Taxonomy tree
# ---------------------------------------------------------------------------

def _usage():
    help(os.path.splitext(os.path.basename(__file__))[0])
    sys.exit(0)


def build_taxonomy_tree(conn):
    """
    Build a parent->children taxonomy tree covering all nodes in dfam_taxdb,
    plus any ancestors needed to connect them to the root.

    Returns {tax_id: {"parent": parent_id | None, "children": [child_ids]}}.
    """
    LOGGER.info("Building taxonomy tree ...")

    rows = list(conn.execute(text(
        "SELECT dfam_taxdb.tax_id, ncbi_taxdb_nodes.parent_id "
        "FROM ncbi_taxdb_nodes "
        "JOIN dfam_taxdb ON dfam_taxdb.tax_id = ncbi_taxdb_nodes.tax_id"
    )))
    tax_ids    = [r[0] for r in rows]
    parent_ids = [r[1] for r in rows]

    while True:
        id_set  = set(tax_ids)
        missing = [p for p in parent_ids if p not in id_set]
        if not missing:
            break
        extra = list(conn.execute(text(
            f"SELECT tax_id, parent_id FROM ncbi_taxdb_nodes "
            f"WHERE tax_id IN ({','.join(str(n) for n in missing)})"
        )))
        for r in extra:
            tax_ids.append(r[0])
            parent_ids.append(r[1])

    T = {tid: {"parent": pid, "children": []} for tid, pid in zip(tax_ids, parent_ids)}
    for n in T:
        p = T[n]["parent"]
        if p in T:
            T[p]["children"].append(n)

    if 1 in T:
        if 1 in T[1]["children"]:
            T[1]["children"].remove(1)
        T[1]["parent"] = None

    LOGGER.info(f"  Tree contains {len(T):,} nodes")
    return T


def _post_order(T, root):
    """Iterative post-order traversal; yields each node after all its descendants."""
    stack = [(root, False)]
    while stack:
        n, visited = stack.pop()
        if visited:
            yield n
        else:
            stack.append((n, True))
            for child in reversed(T[n]["children"]):
                stack.append((child, False))


def subtree_taxa(T, root_id):
    """Return the set of all tax_ids in the subtree rooted at root_id."""
    if root_id not in T:
        LOGGER.warning(f"  taxon {root_id} not found in taxonomy tree")
        return set()
    result = set()
    stack = [root_id]
    while stack:
        n = stack.pop()
        result.add(n)
        stack.extend(T[n]["children"])
    return result


def subtree_taxa_db(conn, taxon_id):
    """
    Return the set of dfam_taxdb tax_ids in the subtree of taxon_id,
    using the dfam_taxdb.lineage field rather than tree traversal.
    """
    row = conn.execute(text(
        "SELECT scientific_name FROM dfam_taxdb WHERE tax_id = :tid"
    ), {"tid": taxon_id}).first()
    if not row:
        row = conn.execute(text(
            "SELECT name_txt FROM ncbi_taxdb_names "
            "WHERE tax_id = :tid AND name_class = 'scientific name'"
        ), {"tid": taxon_id}).first()
    if not row:
        LOGGER.warning(f"  taxon {taxon_id} not found - skipping")
        return set()

    name = row[0]
    rows = conn.execute(text(
        "SELECT tax_id FROM dfam_taxdb "
        "WHERE scientific_name = :name "
        "OR lineage LIKE :mid "
        "OR lineage LIKE :end"
    ), {"name": name, "mid": f"%;{name};%", "end": f"%;{name}"})
    return {int(r[0]) for r in rows}


# ---------------------------------------------------------------------------
# Discovery: top-down greedy cut
# ---------------------------------------------------------------------------

def _per_node_famcounts(conn, accession_prefix):
    """
    Query COUNT(DISTINCT family.id) grouped by dfam_taxdb_tax_id.

    Note: a family with N clade entries is counted once per clade node, so
    summing these over a subtree may overcount.  Acceptable for discovery.
    """
    rows = conn.execute(text(
        "SELECT family_clade.dfam_taxdb_tax_id, COUNT(DISTINCT family.id) "
        "FROM family_clade JOIN family ON family_clade.family_id = family.id "
        f"WHERE family.disabled != 1 AND family.accession LIKE '{accession_prefix}%' "
        "GROUP BY family_clade.dfam_taxdb_tax_id"
    ))
    return {int(r[0]): int(r[1]) for r in rows}


def _build_subtree_counts(T, per_node):
    """
    Bottom-up accumulation of approximate subtree family counts.
    Returns {tax_id: approx_subtree_count}.
    """
    counts = {}
    for n in _post_order(T, 1):
        counts[n] = per_node.get(n, 0) + sum(counts[c] for c in T[n]["children"])
    return counts


def greedy_cut(T, subtree_counts, min_fraction, max_groups):
    """
    Top-down greedy cut: find the most specific subtree roots that individually
    contain >= min_fraction of all families.

    At each node above the threshold, if any child is also above threshold we
    recurse (this node is a waypoint, not a natural boundary).  When no child
    individually clears the bar, the current node is the cut point.

    Returns a list of tax_ids sorted by descending approximate subtree count.
    """
    total = subtree_counts.get(1, 0)
    if total == 0:
        return []
    threshold = total * min_fraction

    groups = []
    stack = list(T[1]["children"])

    while stack:
        n = stack.pop()
        sc = subtree_counts.get(n, 0)
        if sc < threshold:
            continue
        large_children = [c for c in T[n]["children"]
                          if subtree_counts.get(c, 0) >= threshold]
        if large_children:
            stack.extend(T[n]["children"])
        else:
            if len(groups) < max_groups:
                groups.append(n)
            else:
                LOGGER.warning(
                    f"  max-groups limit ({max_groups}) reached; "
                    f"increase --max-groups to find more"
                )
                break

    return sorted(groups, key=lambda x: -subtree_counts.get(x, 0))


def _taxon_names(conn, tax_ids):
    """
    Return {tax_id: {"scientific": str, "common": str|None}}.
    Queries dfam_taxdb first (has curated common names), falls back to
    ncbi_taxdb_names for any IDs not found there.
    """
    if not tax_ids:
        return {}
    id_list = ",".join(str(t) for t in tax_ids)
    names = {}

    rows = conn.execute(text(
        f"SELECT tax_id, scientific_name, common_name FROM dfam_taxdb "
        f"WHERE tax_id IN ({id_list})"
    ))
    for r in rows:
        names[int(r[0])] = {"scientific": r[1], "common": r[2]}

    missing = [t for t in tax_ids if t not in names]
    if missing:
        miss_list = ",".join(str(t) for t in missing)
        rows = conn.execute(text(
            f"SELECT tax_id, name_class, name_txt FROM ncbi_taxdb_names "
            f"WHERE tax_id IN ({miss_list}) "
            f"AND name_class IN ('scientific name', 'common name')"
        ))
        for r in rows:
            tid, cls, txt = int(r[0]), r[1], r[2]
            if tid not in names:
                names[tid] = {"scientific": None, "common": None}
            if cls == "scientific name":
                names[tid]["scientific"] = txt
            elif cls == "common name" and names[tid]["common"] is None:
                names[tid]["common"] = txt

    return names


def _print_suggest_table(label, group_ids, subtree_counts, names):
    total = subtree_counts.get(1, 0)
    bar   = "=" * 74
    print(f"\n{bar}")
    print(f"Suggested {label} groups  "
          f"(min fraction of {total:,} approx annotations)")
    print(bar)
    print(f"  {'Taxon ID':>8}  {'~Count':>12}  {'%':>5}  "
          f"{'Scientific Name':<28}  Common Name")
    print(f"  {'-'*8}  {'-'*12}  {'-'*5}  {'-'*28}  {'-'*28}")
    for tid in group_ids:
        sc  = subtree_counts.get(tid, 0)
        pct = sc / total * 100 if total else 0
        entry = names.get(tid, {})
        sci   = entry.get("scientific") or str(tid)
        cmn   = entry.get("common") or ""
        print(f"  {tid:>8}  {sc:>12,}  {pct:>4.1f}%  {sci:<28}  {cmn}")
    print(bar)


def _print_suggest_snippet(const_name, group_ids, names):
    print(f"\n# ---- copy into famdb_pie_stats.py {const_name} list ----")
    print("# Review/update display names before deploying.")
    print(f"{const_name} = [")
    for tid in group_ids:
        entry  = names.get(tid, {})
        sci    = entry.get("scientific") or str(tid)
        cmn    = entry.get("common")
        # Use common name as display when available; show scientific as comment
        display = cmn if cmn else sci
        comment = f"  # {sci}" if cmn else ""
        print(f'    ("{display:<32}", {tid:>8}),{comment}')
    print("]")
    print(f"# ---- end {const_name} ----\n")


def run_suggest(conn, T, min_fraction, max_groups):
    """
    Run greedy-cut discovery independently for curated and uncurated families.
    Prints tables and Python snippets for both, then returns
    (curated_groups, uncurated_groups) as [(label, tax_id), ...] lists.
    """
    LOGGER.info("Querying per-node family counts for discovery ...")
    curated_pn   = _per_node_famcounts(conn, "DF")
    uncurated_pn = _per_node_famcounts(conn, "DR")

    LOGGER.info("Accumulating subtree counts ...")
    sub_cur   = _build_subtree_counts(T, curated_pn)
    sub_uncur = _build_subtree_counts(T, uncurated_pn)

    for label, sub in [("curated", sub_cur), ("uncurated", sub_uncur)]:
        total = sub.get(1, 0)
        LOGGER.info(
            f"  {label}: approx {total:,} annotations  "
            f"threshold: >= {min_fraction:.1%} ({int(total * min_fraction):,})"
        )

    LOGGER.info("Running greedy cut ...")
    cur_ids   = greedy_cut(T, sub_cur,   min_fraction, max_groups)
    uncur_ids = greedy_cut(T, sub_uncur, min_fraction, max_groups)

    if not cur_ids:
        LOGGER.warning("No curated groups found - try lowering --min-fraction")
    if not uncur_ids:
        LOGGER.error("No uncurated groups found - try lowering --min-fraction")
        sys.exit(1)

    all_ids = list(set(cur_ids) | set(uncur_ids))
    names   = _taxon_names(conn, all_ids)

    _print_suggest_table("CURATED",   cur_ids,   sub_cur,   names)
    _print_suggest_snippet("CURATED_GROUPS",   cur_ids,   names)

    _print_suggest_table("UNCURATED", uncur_ids, sub_uncur, names)
    _print_suggest_snippet("UNCURATED_GROUPS", uncur_ids, names)

    print("NOTE: counts above are approximate (multi-clade families counted once")
    print("      per clade node).  JSON output uses accurate DISTINCT counts.\n")

    def _to_groups(ids, names_dict):
        result = []
        for tid in ids:
            entry   = names_dict.get(tid, {})
            cmn     = entry.get("common")
            sci     = entry.get("scientific") or str(tid)
            result.append((cmn if cmn else sci, tid))
        return result

    return _to_groups(cur_ids, names), _to_groups(uncur_ids, names)


# ---------------------------------------------------------------------------
# Accurate counting (production)
# ---------------------------------------------------------------------------

def count_families_in_group(conn, accession_prefix, taxa_set):
    """
    Count distinct families with any clade in taxa_set.
    For large sets, collects all matching IDs across chunks then deduplicates.
    """
    if not taxa_set:
        return 0
    taxa_list = list(taxa_set)
    if len(taxa_list) <= _CHUNK:
        id_list = ",".join(str(t) for t in taxa_list)
        row = conn.execute(text(
            f"SELECT COUNT(DISTINCT family.id) "
            f"FROM family_clade "
            f"JOIN family ON family_clade.family_id = family.id "
            f"WHERE family.disabled != 1 "
            f"AND family.accession LIKE '{accession_prefix}%' "
            f"AND family_clade.dfam_taxdb_tax_id IN ({id_list})"
        )).one()
        return int(row[0])

    matching_ids = set()
    for i in range(0, len(taxa_list), _CHUNK):
        chunk   = taxa_list[i:i + _CHUNK]
        id_list = ",".join(str(t) for t in chunk)
        rows = conn.execute(text(
            f"SELECT DISTINCT family.id "
            f"FROM family_clade "
            f"JOIN family ON family_clade.family_id = family.id "
            f"WHERE family.disabled != 1 "
            f"AND family.accession LIKE '{accession_prefix}%' "
            f"AND family_clade.dfam_taxdb_tax_id IN ({id_list})"
        ))
        for row in rows:
            matching_ids.add(row[0])
    return len(matching_ids)


def count_other_families(conn, accession_prefix, all_named_taxa):
    """Count families with no clade in any named group (total minus those in any group)."""
    (total,) = conn.execute(text(
        f"SELECT COUNT(DISTINCT id) FROM family "
        f"WHERE disabled != 1 AND accession LIKE '{accession_prefix}%'"
    )).one()
    total = int(total)
    if not all_named_taxa:
        return total
    return total - count_families_in_group(conn, accession_prefix, all_named_taxa)


_SPECIES_RANKS = (
    "'species', 'subspecies', 'varietas', 'forma', 'forma specialis', 'no rank', 'strain'"
)

def count_species_in_group(conn, accession_prefix, taxa_set):
    """
    Count distinct species-level taxa with at least one family clade within taxa_set.
    Counts taxa whose ncbi_taxdb_nodes.rank is in the species-level rank set.
    """
    if not taxa_set:
        return 0
    taxa_list = list(taxa_set)
    if len(taxa_list) <= _CHUNK:
        id_list = ",".join(str(t) for t in taxa_list)
        row = conn.execute(text(
            f"SELECT COUNT(DISTINCT family_clade.dfam_taxdb_tax_id) "
            f"FROM family_clade "
            f"JOIN family ON family_clade.family_id = family.id "
            f"JOIN ncbi_taxdb_nodes ON family_clade.dfam_taxdb_tax_id = ncbi_taxdb_nodes.tax_id "
            f"WHERE family.disabled != 1 "
            f"AND family.accession LIKE '{accession_prefix}%' "
            f"AND ncbi_taxdb_nodes.rank IN ({_SPECIES_RANKS}) "
            f"AND family_clade.dfam_taxdb_tax_id IN ({id_list})"
        )).one()
        return int(row[0])

    matching_taxa = set()
    for i in range(0, len(taxa_list), _CHUNK):
        chunk   = taxa_list[i:i + _CHUNK]
        id_list = ",".join(str(t) for t in chunk)
        rows = conn.execute(text(
            f"SELECT DISTINCT family_clade.dfam_taxdb_tax_id "
            f"FROM family_clade "
            f"JOIN family ON family_clade.family_id = family.id "
            f"JOIN ncbi_taxdb_nodes ON family_clade.dfam_taxdb_tax_id = ncbi_taxdb_nodes.tax_id "
            f"WHERE family.disabled != 1 "
            f"AND family.accession LIKE '{accession_prefix}%' "
            f"AND ncbi_taxdb_nodes.rank IN ({_SPECIES_RANKS}) "
            f"AND family_clade.dfam_taxdb_tax_id IN ({id_list})"
        ))
        for row in rows:
            matching_taxa.add(row[0])
    return len(matching_taxa)


def count_other_species(conn, accession_prefix, all_named_taxa):
    """Count species-level taxa with families but no clade in any named group."""
    (total,) = conn.execute(text(
        f"SELECT COUNT(DISTINCT family_clade.dfam_taxdb_tax_id) "
        f"FROM family_clade "
        f"JOIN family ON family_clade.family_id = family.id "
        f"JOIN ncbi_taxdb_nodes ON family_clade.dfam_taxdb_tax_id = ncbi_taxdb_nodes.tax_id "
        f"WHERE family.disabled != 1 "
        f"AND family.accession LIKE '{accession_prefix}%' "
        f"AND ncbi_taxdb_nodes.rank IN ({_SPECIES_RANKS})"
    )).one()
    total = int(total)
    if not all_named_taxa:
        return total
    return total - count_species_in_group(conn, accession_prefix, all_named_taxa)


def _log_other_candidates(conn, T, accession_prefix, all_named_taxa, top_n=5):
    """Log the largest approximate-count candidate subtrees within the Other bucket."""
    per_node = _per_node_famcounts(conn, accession_prefix)

    # Build subtree counts treating named-group nodes as walls (zeroed out)
    counts = {}
    for n in _post_order(T, 1):
        if n in all_named_taxa:
            counts[n] = 0
        else:
            counts[n] = per_node.get(n, 0) + sum(counts.get(c, 0) for c in T[n]["children"])

    total_other_approx = counts.get(1, 0)
    if total_other_approx == 0:
        return

    # Top-down greedy: descend while children account for most of the parent's count
    candidates = []
    stack = [1]
    while stack:
        n = stack.pop()
        if n in all_named_taxa:
            stack.extend(T[n]["children"])
            continue
        nc = counts.get(n, 0)
        if nc == 0:
            continue
        child_sum = sum(counts.get(c, 0) for c in T[n]["children"])
        if child_sum >= nc * 0.8:
            stack.extend(T[n]["children"])
        else:
            candidates.append((nc, n))

    candidates.sort(reverse=True)
    top = candidates[:top_n]
    if not top:
        return

    names = _taxon_names(conn, [n for _, n in top])
    LOGGER.info(f"  Largest sub-groups within Other (approximate counts):")
    for approx_count, tid in top:
        entry = names.get(tid, {})
        name = entry.get("common") or entry.get("scientific") or str(tid)
        pct = approx_count / total_other_approx * 100
        LOGGER.info(f"    {name} (taxon {tid}): ~{approx_count:,}  [{pct:.0f}% of Other]")


def compute_stats(conn, T, accession_prefix, label, groups, count_by="families"):
    """
    For one curation class (DF* or DR*), compute per-group accurate counts.
    count_by is "families" (distinct TE families) or "species" (distinct species-rank taxa).
    Returns a list of {"group", "taxon", "count"} dicts; zero-count groups omitted;
    "Other" appended last if non-zero.
    """
    LOGGER.info(f"Computing {label} {count_by} counts ...")
    all_named_taxa = set()
    results = []

    for group_name, taxon_id in groups:
        taxa = subtree_taxa_db(conn, taxon_id)
        all_named_taxa |= taxa
        if count_by == "species":
            count = count_species_in_group(conn, accession_prefix, taxa)
        else:
            count = count_families_in_group(conn, accession_prefix, taxa)
        LOGGER.info(f"  {group_name}: {count:,}")
        if count > 0:
            results.append({"group": group_name, "taxon": taxon_id, "count": count})

    if count_by == "species":
        other = count_other_species(conn, accession_prefix, all_named_taxa)
    else:
        other = count_other_families(conn, accession_prefix, all_named_taxa)
    LOGGER.info(f"  Other: {other:,}")
    if other > 0:
        results.append({"group": "Other", "taxon": 0, "count": other})
        if count_by == "families":
            _log_other_candidates(conn, T, accession_prefix, all_named_taxa)

    return results


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    logging.basicConfig(stream=sys.stdout, format="%(levelname)s: %(message)s")

    class _HelpAction(argparse.Action):
        def __init__(self, option_strings, dest, default=False, required=False, help=None):
            super().__init__(option_strings=option_strings, dest=dest, nargs=0,
                             const=True, default=default, required=required, help=help)
        def __call__(self, parser, args, values, option_string=None):
            _usage()

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action=_HelpAction)
    parser.add_argument("-l", "--log-level", default="INFO")
    parser.add_argument("-c", "--dfam-config")

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    subparsers.required = True

    gen_parser = subparsers.add_parser(
        "generate",
        help="Compute per-group family counts and write the pie-stats JSON",
        add_help=True,
    )
    gen_parser.add_argument("-g", "--groups-file", metavar="JSON",
                            help="JSON file with curated/uncurated taxon lists "
                                 "(same format as this script's output); counts are "
                                 "ignored and recalculated")
    gen_parser.add_argument("-o", "--output", default="pie_stats.json",
                            help="Output JSON file (default: pie_stats.json)")
    gen_parser.add_argument("-b", "--count-by", choices=["families", "species"],
                            default="families",
                            help="Count distinct families (default) or species-rank taxa")

    sug_parser = subparsers.add_parser(
        "suggest",
        help="Discover candidate group lists from the live database",
        add_help=True,
    )
    sug_parser.add_argument("-o", "--output", default="pie_stats.json",
                            help="Output JSON file (default: pie_stats.json)")
    sug_parser.add_argument("-f", "--min-fraction", type=float, default=0.03, metavar="F",
                            help="Min fraction of families for a group (default: 0.03)")
    sug_parser.add_argument("-m", "--max-groups", type=int, default=20, metavar="N",
                            help="Max groups to suggest per dataset (default: 20)")
    sug_parser.add_argument("-b", "--count-by", choices=["families", "species"],
                            default="families",
                            help="Count distinct families (default) or species-rank taxa "
                                 "in the final JSON output (discovery always uses families)")

    args = parser.parse_args()

    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

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

    conf   = dc.DfamConfig(args.dfam_config)
    dfamdb = create_engine(conf.getDBConnStrWPassFallback("Dfam"))
    session = sessionmaker(dfamdb)()

    version_info = session.execute(select(dfam.DbVersion)).scalar_one()
    db_version = version_info.dfam_version
    db_date    = version_info.dfam_release_date.strftime("%Y-%m-%d")
    LOGGER.info(f"Dfam {db_version} ({db_date})")

    db_url = dfamdb.url
    LOGGER.info(
        f"Connected to '{db_url.database}' on {db_url.host}"
        + (f":{db_url.port}" if db_url.port else "")
        + f" as '{db_url.username}'"
    )

    with dfamdb.connect() as conn:
        T = build_taxonomy_tree(conn)

        if args.command == "suggest":
            cur_groups, uncur_groups = run_suggest(
                conn, T, args.min_fraction, args.max_groups
            )
            LOGGER.info("Computing accurate counts for suggested groups ...")
        else:
            if args.groups_file:
                with open(args.groups_file) as gf:
                    gdata = json.load(gf)
                def _load_groups(entries):
                    return [
                        (e["group"], e["taxon"])
                        for e in entries
                        if e.get("taxon", 0) != 0
                    ]
                cur_groups   = _load_groups(gdata.get("curated",   []))
                uncur_groups = _load_groups(gdata.get("uncurated", []))
                LOGGER.info(
                    f"Loaded {len(cur_groups)} curated and {len(uncur_groups)} "
                    f"uncurated groups from {args.groups_file}"
                )
            else:
                cur_groups   = CURATED_GROUPS
                uncur_groups = UNCURATED_GROUPS

        count_by = args.count_by
        curated   = compute_stats(conn, T, "DF", "curated (DF*)",   cur_groups,   count_by)
        uncurated = compute_stats(conn, T, "DR", "uncurated (DR*)", uncur_groups, count_by)

    output = {
        "curated":   curated,
        "uncurated": uncurated,
        "meta": {
            "db_version": db_version,
            "db_date":    db_date,
            "count_by":   count_by,
        },
    }

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)
    LOGGER.info(f"Written to {args.output}")


if __name__ == "__main__":
    main()
