"""
Fakes, stubs, etc. for use in testing FamDB v3
"""

from copy import deepcopy
from famdb_classes import FamDBLeaf, FamDBRoot
from famdb_helper_classes import TaxNode, Family
from famdb_globals import (
    FAMDB_VERSION,
    DESCRIPTION,
    META_META,
    META_UUID,
    META_DB_VERSION,
    META_DB_DATE,
    META_FILE_MAP,
    META_FAMDB_VERSION,
    META_CREATED,
    META_DB_DESCRIPTION,
    COMPONENT_CC,
    COMPONENT_CH,
    COMPONENT_UC,
    COMPONENT_UH,
)

"""
Taxonomy tree used in tests:

        1
      /   \\
 (0) 2     3
--------------
(1)/ |\\ (2)
  4  |  *5
 /   |   \\
6    |    7
     (2)

Partition assignments (node -> component partitions):
  NODES_CC = {0: [1,2,3,4,5,6,7]}  (all taxa, single curated-consensus file)
  NODES_CH = {1: [1,2,3], 2: [4,5,6,7]}
  NODES_UC = {1: [5,6,7]}
"""

TAX_NAMES = {
    1: "root",
    2: "Order",
    3: "Other Order",
    4: "Genus",
    5: "Other Genus",
    6: "Species",
    7: "Other Species",
}
COMMON_NAMES = {
    1: "Root Dummy 1",
    2: "Root Dummy 2",
    3: "Root Dummy 3",
    4: "Leaf Dummy 4",
    5: "Leaf Dummy 5",
    6: "Leaf Dummy 6",
    7: "Leaf Dummy 7",
}

# Old-style partition nodes (kept for init_single_file compatibility)
NODES = {0: [1, 2, 3], 1: [4, 6], 2: [5, 7]}

# Component-aware partition node lists
NODES_CC = {0: [1, 2, 3, 4, 5, 6, 7]}   # all taxa in single CC file
NODES_CH = {1: [1, 2, 3], 2: [4, 5, 6, 7]}
NODES_UC = {1: [5, 6, 7]}

FILE_INFO = {
    META_META: {META_UUID: "uuidXX", META_DB_VERSION: "V1", META_DB_DATE: "2020-07-15"},
    META_FILE_MAP: {
        "0": {
            "T_root": 1,
            "filename": "unittest.0.h5",
            "F_roots": [],
            "T_root_name": "Root Node",
            "F_roots_names": [],
        },
        "cc.0": {
            "T_root": 1,
            "filename": "unittest.curated.consensus.0.h5",
            "F_roots": [1],
            "T_root_name": "Root Node",
            "F_roots_names": [],
        },
        "ch.1": {
            "T_root": 1,
            "filename": "unittest.curated.hmm.1.h5",
            "F_roots": [1],
            "T_root_name": "Root Node",
            "F_roots_names": [],
        },
        "ch.2": {
            "T_root": 4,
            "filename": "unittest.curated.hmm.2.h5",
            "F_roots": [4],
            "T_root_name": "Genus Node",
            "F_roots_names": [],
        },
        "uc.1": {
            "T_root": 5,
            "filename": "unittest.uncurated.consensus.1.h5",
            "F_roots": [5],
            "T_root_name": "Other Genus",
            "F_roots_names": [],
        },
    },
}

DB_INFO = ("Test Dfam", "V1", "2020-07-15", "<copyright header>")
FAKE_REPPEPS = "./tests/rep_pep_test.lib"


def build_taxa(nodes):
    for node in nodes.values():
        if node.tax_id != 1:
            node.parent_node = nodes[node.parent_id]
            node.parent_node.children += [node]
        node.names += [["scientific name", TAX_NAMES[node.tax_id]]]
        node.names += [["common name", COMMON_NAMES[node.tax_id]]]
    return nodes


def make_family(acc, clades, consensus, model):
    """Convenience factory to generate a test Family object."""
    fam = Family()
    fam.accession = acc
    fam.name = "Test family " + acc
    fam.version = 1
    fam.clades = clades
    fam.consensus = consensus
    fam.model = model
    return fam


def write_test_metadata(db):
    """Override format metadata for testing (bypasses __write_metadata timestamp)."""
    db.file.attrs[META_FAMDB_VERSION] = FAMDB_VERSION
    db.file.attrs[META_CREATED] = "<creation date>"
    db.file.attrs[META_DB_DESCRIPTION] = DESCRIPTION


def _make_partition_cache(taxa):
    """Build the PartitionCache dict for the test taxonomy."""
    cache = {}
    # Curated families: TEST000x are in taxa 1-4 (CC.0, CH varies)
    # Uncurated families: DR* are in taxa 6,7 (UC.1)
    curated_taxa = {1, 2, 3, 4}
    uncurated_taxa = {5, 6, 7}
    for tax_id in taxa:
        entry = {
            COMPONENT_CC: None,
            COMPONENT_CH: None,
            COMPONENT_UC: None,
            COMPONENT_UH: None,
        }
        if tax_id in curated_taxa:
            entry[COMPONENT_CC] = 0
            entry[COMPONENT_CH] = 1 if tax_id in {1, 2, 3} else 2
        if tax_id in uncurated_taxa:
            entry[COMPONENT_UC] = 1
        # Only store taxa that have at least one component
        if any(v is not None for v in entry.values()):
            cache[str(tax_id)] = entry
    return cache


def _make_family_taxon_map(families):
    """Build the Lookup/ByTaxon mapping from a list of Family objects."""
    mapping = {}
    for fam in families:
        for clade_id in fam.clades:
            mapping.setdefault(str(clade_id), []).append(fam.accession)
    return mapping


def init_db_file(filename):
    """
    Creates a v3 FamDB test database with the following files:
      {filename}.0.h5                       — root
      {filename}.curated.consensus.0.h5     — curated consensus (CC)
      {filename}.curated.hmm.1.h5           — curated HMM partition 1
      {filename}.curated.hmm.2.h5           — curated HMM partition 2
      {filename}.uncurated.consensus.1.h5   — uncurated consensus (UC)
    """
    # Curated families: TEST* prefix (treated as curated for testing)
    CURATED = [
        make_family("TEST0001", [1], "ACGT", "<model1>"),
        make_family("TEST0002", [2, 3], None, "<model2>"),
        make_family("TEST0003", [3], "GGTC", "<model3>"),
        make_family("TEST0004", [4], "CCCCTTTT", None),
    ]
    # Uncurated families: DR* prefix
    UNCURATED = [
        make_family("DR000000001", [7], "GCATATCG", None),
        make_family("DR_Repeat1", [6], "CGACTAT", None),
    ]

    CURATED[1].name = None          # TEST0002 has no name
    CURATED[2].search_stages = "30,40"
    CURATED[3].buffer_stages = "10[1-2],10[5-8],20"
    CURATED[3].search_stages = "35"
    CURATED[3].repeat_type = "SINE"

    ALL_FAMILIES = CURATED + UNCURATED

    TAX_DB = {
        1: TaxNode(1, None),
        2: TaxNode(2, 1),
        3: TaxNode(3, 1),
        4: TaxNode(4, 2),
        5: TaxNode(5, 2),
        6: TaxNode(6, 4),
        7: TaxNode(7, 5),
    }
    taxa = build_taxa(TAX_DB)

    partition_cache = _make_partition_cache(TAX_DB.keys())
    family_taxon_map = _make_family_taxon_map(ALL_FAMILIES)

    # --- Root file ---
    with FamDBRoot(f"{filename}.0.h5", "w") as db:
        db.set_metadata("0", FILE_INFO, *DB_INFO, is_root=True)
        write_test_metadata(db)
        db.write_repeatpeps(FAKE_REPPEPS)
        db.write_full_taxonomy(taxa)
        db.write_lookup_bytaxon(family_taxon_map)
        db.write_partition_cache(partition_cache)
        db.finalize()

    # --- Curated consensus (CC.0): all curated families ---
    with FamDBLeaf(f"{filename}.curated.consensus.0.h5", "w", component_type=COMPONENT_CC) as db:
        db.set_metadata("cc.0", FILE_INFO, *DB_INFO)
        write_test_metadata(db)
        for fam in CURATED:
            db.add_family(fam)
        db.finalize()

    # --- Curated HMM partition 1: taxa [1,2,3] ---
    ch1_fams = [f for f in CURATED if any(c in {1, 2, 3} for c in f.clades)]
    with FamDBLeaf(f"{filename}.curated.hmm.1.h5", "w", component_type=COMPONENT_CH) as db:
        db.set_metadata("ch.1", FILE_INFO, *DB_INFO)
        write_test_metadata(db)
        for fam in ch1_fams:
            db.add_family(fam)
        db.finalize()

    # --- Curated HMM partition 2: taxa [4,5,6,7] ---
    ch2_fams = [f for f in CURATED if any(c in {4, 5, 6, 7} for c in f.clades) and
                not any(c in {1, 2, 3} for c in f.clades)]
    with FamDBLeaf(f"{filename}.curated.hmm.2.h5", "w", component_type=COMPONENT_CH) as db:
        db.set_metadata("ch.2", FILE_INFO, *DB_INFO)
        write_test_metadata(db)
        for fam in ch2_fams:
            db.add_family(fam)
        db.finalize()

    # --- Uncurated consensus partition 1: all uncurated families ---
    with FamDBLeaf(f"{filename}.uncurated.consensus.1.h5", "w", component_type=COMPONENT_UC) as db:
        db.set_metadata("uc.1", FILE_INFO, *DB_INFO)
        write_test_metadata(db)
        for fam in UNCURATED:
            db.add_family(fam)
        db.finalize()


def init_single_file(n, db_dir, change_id=False):
    """
    Creates a single file for validation tests.
    Only creates the root file (n=0) or a legacy-style leaf file.
    Used by file_checker tests.
    """
    TAX_DB = {
        1: TaxNode(1, None),
        2: TaxNode(2, 1),
        3: TaxNode(3, 1),
        4: TaxNode(4, 2),
        5: TaxNode(5, 2),
        6: TaxNode(6, 4),
        7: TaxNode(7, 5),
    }
    taxa = build_taxa(TAX_DB)

    if change_id:
        file_info = deepcopy(FILE_INFO)
        file_info[META_META][META_UUID] = "uuidYY"
    else:
        file_info = deepcopy(FILE_INFO)

    if n == 0:
        filename = f"{db_dir}.0.h5"
        partition_cache = _make_partition_cache(TAX_DB.keys())
        with FamDBRoot(filename, "w") as f:
            f.write_full_taxonomy(taxa)
            f.write_lookup_bytaxon({})
            f.write_partition_cache(partition_cache)
            write_test_metadata(f)
            f.set_metadata("0", file_info, *DB_INFO, is_root=True)
            f.finalize()
    else:
        # Create a minimal curated consensus file for validation
        filename = f"{db_dir}.curated.consensus.{n}.h5"
        with FamDBLeaf(filename, "w", component_type=COMPONENT_CC) as f:
            write_test_metadata(f)
            f.set_metadata(f"cc.{n}", file_info, *DB_INFO)
            f.finalize()
