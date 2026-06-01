import os
import unittest
from famdb_classes import FamDBLeaf, FamDBRoot, FamDB
from famdb_helper_classes import Family
from .doubles import init_db_file, FILE_INFO
from unittest.mock import patch
import io
from famdb_globals import FAMDB_VERSION, TEST_DIR, DESCRIPTION, COMPONENT_CC, COMPONENT_CH, COMPONENT_UC, COMPONENT_UH


# Partition cache dicts for each taxon group in the test fixture:
#   Taxa 1, 2, 3 → curated, ch.1
#   Taxon 4      → curated, ch.2
#   Taxa 5, 6, 7 → uncurated, uc.1
PC_HIGH = {"cc": 0, "ch": 1, "uc": None, "uh": None}  # taxa 1, 2, 3
PC_LOW  = {"cc": 0, "ch": 2, "uc": None, "uh": None}  # taxon 4
PC_UC   = {"cc": None, "ch": None, "uc": 1, "uh": None}  # taxa 5, 6, 7


class TestDatabase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        file_dir = f"{TEST_DIR}/db"
        os.makedirs(file_dir, exist_ok=True)
        db_dir = f"{file_dir}/unittest"
        init_db_file(db_dir)
        filenames = [
            f"{db_dir}.0.h5",
            f"{db_dir}.curated.consensus.0.h5",
            f"{db_dir}.curated.hmm.1.h5",
            f"{db_dir}.curated.hmm.2.h5",
            f"{db_dir}.uncurated.consensus.1.h5",
        ]
        TestDatabase.filenames = filenames
        TestDatabase.file_dir = file_dir
        TestDatabase.famdb = FamDB(file_dir, "r+")
        TestDatabase.famdb.build_pruned_tree()

    @classmethod
    def tearDownClass(cls):
        filenames = TestDatabase.filenames
        TestDatabase.filenames = None

        for name in filenames:
            if os.path.exists(name):
                os.remove(name)
        os.rmdir(TestDatabase.file_dir)

    def test_get_metadata(self):
        # CC.0 leaf: partition "cc.0" → T_root_name "Root Node", F_roots_names []
        test_info = {
            "famdb_version": FAMDB_VERSION,
            "created": "<creation date>",
            "partition_name": "Root Node",
            "partition_detail": "",
            "name": "Test Dfam",
            "db_version": "V1",
            "date": "2020-07-15",
            "description": DESCRIPTION,
            "copyright": "<copyright header>",
        }

        with FamDBLeaf(TestDatabase.filenames[1], "r") as db:
            self.assertEqual(db.get_metadata(), test_info)

        # Root file: partition "0" → T_root_name "Root Node"
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(
                db.get_metadata(),
                test_info,
            )

    def test_get_history(self):
        substrings = [
            "File Initialized",
            "Metadata Set",
            "RepeatPeps Written",
            "Taxonomy Nodes Written",
        ]
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            history = db.get_history()
            for substring in substrings:
                self.assertIn(substring, history)

    def test_interrupt_check(self):
        message = "Test Message"
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            stamp = db.update_changelog(message)
            self.assertTrue(db.interrupt_check())
            db._verify_change(stamp, message)
            self.assertFalse(db.interrupt_check())

    def test_update_description(self):
        new_desc = "New Description"
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            db.update_description(new_desc)
            self.assertEqual(db.get_metadata()["description"], new_desc)

    def test_get_counts(self):
        # Root has no families
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_counts(), {"consensus": 0, "hmm": 0})

        # CC.0: TEST0001(cons), TEST0002(no cons), TEST0003(cons), TEST0004(cons) → 3 consensus
        with FamDBLeaf(TestDatabase.filenames[1], "r") as db:
            self.assertEqual(db.get_counts(), {"consensus": 3, "hmm": 0})

        # CH.1: TEST0001(model), TEST0002(model), TEST0003(model) → 3 hmm
        with FamDBLeaf(TestDatabase.filenames[2], "r") as db:
            self.assertEqual(db.get_counts(), {"consensus": 0, "hmm": 3})

        # CH.2: TEST0004(no model) → 0
        with FamDBLeaf(TestDatabase.filenames[3], "r") as db:
            self.assertEqual(db.get_counts(), {"consensus": 0, "hmm": 0})

        # UC.1: DR000000001(cons), DR_Repeat1(cons) → 2 consensus
        with FamDBLeaf(TestDatabase.filenames[4], "r") as db:
            self.assertEqual(db.get_counts(), {"consensus": 2, "hmm": 0})

    def test_get_partition_num(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_partition_num(), "0")

        with FamDBLeaf(TestDatabase.filenames[1], "r") as db:
            self.assertEqual(db.get_partition_num(), "cc.0")

        with FamDBLeaf(TestDatabase.filenames[2], "r") as db:
            self.assertEqual(db.get_partition_num(), "ch.1")

    def test_get_file_info(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertDictEqual(db.get_file_info(), FILE_INFO)

    def test_is_root(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.is_root(), True)

        with FamDBLeaf(TestDatabase.filenames[1], "r") as db:
            self.assertEqual(db.is_root(), False)

    def test_get_family_by_accession(self):
        # Families live in component files, not root
        with FamDBLeaf(TestDatabase.filenames[1], "r") as db:  # CC.0
            test_fam = db.get_family_by_accession("TEST0001")
            self.assertIsInstance(test_fam, Family)
            self.assertEqual(test_fam.name, "Test family TEST0001")
            self.assertEqual(db.get_family_by_accession("TEST0000"), None)

    def test_get_family_by_name(self):
        # Root has no families → name lookup returns None
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_family_by_name("Test family TEST0002"), None)
        # CC.0 has TEST0004
        with FamDBLeaf(TestDatabase.filenames[1], "r") as db:
            test_fam = db.get_family_by_name("Test family TEST0004")
            self.assertIsInstance(test_fam, Family)
            self.assertEqual(test_fam.name, "Test family TEST0004")

    def test_get_families_for_taxon(self):
        # Lookup/ByTaxon lives only in the root file
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_families_for_taxon(3), ["TEST0002", "TEST0003"])
            self.assertEqual(db.get_families_for_taxon(4), ["TEST0004"])

    # Root File Methods ------------------------------------------------
    def test_get_complete_lineage(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_lineage(4), [4])
            self.assertEqual(
                db.get_lineage(4, descendants=True, complete=True), [4, [6]]
            )
            self.assertEqual(
                db.get_lineage(6, ancestors=True, complete=True),
                [1, [2, [4, [6]]]],
            )
            self.assertEqual(
                db.get_lineage(4, ancestors=True, descendants=True, complete=True),
                [1, [2, [4, [6]]]],
            )

            self.assertEqual(db.get_lineage(1, complete=True), [1])
            self.assertEqual(
                db.get_lineage(1, descendants=True, complete=True),
                [1, [2, [4, [6]], [5, [7]]], [3]],
            )
            self.assertEqual(
                db.get_lineage(2, ancestors=True, descendants=True, complete=True),
                [1, [2, [4, [6]], [5, [7]]]],
            )

            self.assertEqual(
                db.get_lineage(5, descendants=True, complete=False), [5, [7]]
            )
            self.assertEqual(
                db.get_lineage(7, ancestors=True, complete=False),
                [1, [2, [7]]],
            )
            self.assertEqual(
                db.get_lineage(5, ancestors=True, complete=False),
                [1, [2, [5]]],
            )

            self.assertEqual(
                db.get_lineage(1, descendants=True, complete=False),
                [1, [2, [4, [6]], [7]], [3]],
            )
            self.assertEqual(
                db.get_lineage(3, ancestors=True, complete=False), [1, [3]]
            )
            self.assertEqual(
                db.get_lineage(2, ancestors=True, descendants=True, complete=False),
                [1, [2, [4, [6]], [7]]],
            )

    def test_search_taxon_names(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(
                list(db.search_taxon_names("Order")),
                [
                    [2, True, PC_HIGH],
                    [3, False, PC_HIGH],
                ],
            )

            self.assertEqual(
                list(db.search_taxon_names("Genus")),
                [
                    [4, True, PC_LOW],
                    [5, False, PC_UC],
                ],
            )

            self.assertEqual(
                list(db.search_taxon_names("rut", search_similar=True)),
                [
                    [1, False, PC_HIGH],
                ],
            )

            self.assertEqual(
                list(db.search_taxon_names("Root Dummy", "common name")),
                [
                    [1, False, PC_HIGH],
                    [2, False, PC_HIGH],
                    [3, False, PC_HIGH],
                ],
            )

            self.assertEqual(list(db.search_taxon_names("Missing")), [])

    def test_get_taxon_name(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_taxon_name(2), ["Order", PC_HIGH])
            self.assertEqual(db.get_taxon_name(10), ("Not Found", "N/A"))
            self.assertEqual(db.get_taxon_name(2, "common name"), ["Root Dummy 2", PC_HIGH])
            self.assertEqual(db.get_taxon_name(4), ["Genus", PC_LOW])

    def test_get_taxon_names(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(
                db.get_taxon_names(2),
                [["scientific name", "Order"], ["common name", "Root Dummy 2"]],
            )
            self.assertEqual(
                db.get_taxon_names(4),
                [["scientific name", "Genus"], ["common name", "Leaf Dummy 4"]],
            )
            self.assertEqual(db.get_taxon_names(10), [])

    def test_get_lineage_path(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(
                db.get_lineage_path(3, complete=True),
                [["root", PC_HIGH], ["Other Order", PC_HIGH]],
            )

            # test caching in get_lineage_path
            self.assertEqual(
                db.get_lineage_path(3, complete=True),
                [["root", PC_HIGH], ["Other Order", PC_HIGH]],
            )

            # test lookup without cache
            self.assertEqual(
                db.get_lineage_path(3, cache=False, complete=True),
                [["root", PC_HIGH], ["Other Order", PC_HIGH]],
            )

            # test with supplied tree
            self.assertEqual(
                db.get_lineage_path(4, [1, [2, [4]], [3]], complete=True),
                [["root", PC_HIGH], ["Order", PC_HIGH], ["Genus", PC_LOW]],
            )

    def test_resolve_species(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.resolve_species(3), [[3, PC_HIGH, True]])
            self.assertEqual(db.resolve_species(4), [[4, PC_LOW, True]])
            self.assertEqual(db.resolve_species(999), [])
            self.assertEqual(
                db.resolve_species("Species"),
                [[6, PC_UC, True], [7, PC_UC, False]],
            )
            self.assertEqual(db.resolve_species("Tardigrade"), [])

    def test_resolve_one_species(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.resolve_one_species(3), [3, PC_HIGH])
            self.assertEqual(db.resolve_one_species(999), (None, None))
            self.assertEqual(db.resolve_one_species("Species"), [6, PC_UC])
            self.assertEqual(db.resolve_one_species("Mus musculus"), (None, None))

    def test_get_sanitized_name(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_sanitized_name(5), "Other_Genus")

    def test_find_files(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_file_info(), FILE_INFO)

    def test_find_taxon(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.find_taxon(2), PC_HIGH)
            self.assertEqual(db.find_taxon(4), PC_LOW)
            self.assertEqual(db.find_taxon(5), PC_UC)

    def test_repeatpeps(self):
        with FamDBRoot(TestDatabase.filenames[0], "r") as db:
            self.assertEqual(db.get_repeatpeps(), ">DUMMYACC\nDUMMYDATA")

    # Umbrella Methods -----------------------------------------------------------------------------
    def test_get_lineage(self):
        famdb = TestDatabase.famdb
        self.assertEqual(famdb.get_lineage(4), [4])
        self.assertEqual(
            famdb.get_lineage(4, descendants=True, complete=True), [4, [6]]
        )
        self.assertEqual(
            famdb.get_lineage(6, ancestors=True, complete=True),
            [1, [2, [4, [6]]]],
        )
        self.assertEqual(
            famdb.get_lineage(4, ancestors=True, descendants=True, complete=True),
            [1, [2, [4, [6]]]],
        )

        self.assertEqual(famdb.get_lineage(1, complete=True), [1])
        self.assertEqual(
            famdb.get_lineage(1, descendants=True, complete=True),
            [1, [2, [4, [6]], [5, [7]]], [3]],
        )
        self.assertEqual(
            famdb.get_lineage(2, ancestors=True, descendants=True, complete=True),
            [1, [2, [4, [6]], [5, [7]]]],
        )

        self.assertEqual(
            famdb.get_lineage(5, descendants=True, complete=False), [5, [7]]
        )
        self.assertEqual(
            famdb.get_lineage(7, ancestors=True, complete=False),
            [1, [2, [7]]],
        )
        self.assertEqual(
            famdb.get_lineage(5, ancestors=True, complete=False),
            [1, [2, [5]]],
        )

        self.assertEqual(
            famdb.get_lineage(1, descendants=True, complete=False),
            [1, [2, [4, [6]], [7]], [3]],
        )
        self.assertEqual(famdb.get_lineage(3, ancestors=True, complete=False), [1, [3]])
        self.assertEqual(
            famdb.get_lineage(2, ancestors=True, descendants=True, complete=False),
            [1, [2, [4, [6]], [7]]],
        )

    @patch("sys.stdout", new_callable=io.StringIO)
    def test_show_files(self, mock_print):
        famdb = TestDatabase.famdb
        famdb.show_files()
        out = (
            "\nInstalled Components\n--------------------\n"
            "\n Curated Consensus:\n"
            "     partition 0 [unittest.curated.consensus.0.h5]:  Root Node  3 families\n"
            "\n Curated HMMs:\n"
            "     partition 1 [unittest.curated.hmm.1.h5]:  Root Node   3 families\n"
            "     partition 2 [unittest.curated.hmm.2.h5]:  Genus Node  0 families\n"
            "\n Uncurated Consensus:\n"
            "     partition 1 [unittest.uncurated.consensus.1.h5]:  Other Genus  2 families\n"
            "\n Uncurated HMMs:\n"
            "     [ Not Installed ]\n"
            "\n"
        )
        self.assertEqual(mock_print.getvalue(), out)

    def test_get_complete_lineage_path(self):
        famdb = TestDatabase.famdb
        self.assertEqual(
            famdb.get_lineage_path(5, cache=False, complete=True),
            [["root", PC_HIGH], ["Order", PC_HIGH], ["Other Genus", PC_UC]],
        )
        self.assertEqual(
            famdb.get_lineage_path(5, partition=False, cache=False, complete=True),
            ["root", "Order", "Other Genus"],
        )

    def test_get_pruned_lineage_path(self):
        famdb = TestDatabase.famdb
        self.assertEqual(
            famdb.get_lineage_path(7, complete=False),
            [["root", PC_HIGH], ["Order", PC_HIGH], ["Other Species", PC_UC]],
        )
        self.assertEqual(
            famdb.get_lineage_path(5, complete=False),
            [["root", PC_HIGH], ["Order", PC_HIGH], ["Other Genus", PC_UC]],
        )
        self.assertEqual(
            famdb.get_lineage_path(5, partition=False, cache=False, complete=False),
            ["root", "Order", "Other Genus"],
        )

    def test_get_counts(self):
        famdb = TestDatabase.famdb
        # root(0,0) + cc.0(3,0) + ch.1(0,3) + ch.2(0,0) + uc.1(2,0) = {c:5, h:3, file:5}
        self.assertEqual(famdb.get_counts(), {"consensus": 5, "hmm": 3, "file": 5})

    def test_resolve_names(self):
        famdb = TestDatabase.famdb
        self.assertEqual(
            famdb.resolve_names(4),
            [
                [
                    4,
                    True,
                    PC_LOW,
                    [["scientific name", "Genus"], ["common name", "Leaf Dummy 4"]],
                ]
            ],
        )
        self.assertEqual(
            famdb.resolve_names(2),
            [
                [
                    2,
                    True,
                    PC_HIGH,
                    [["scientific name", "Order"], ["common name", "Root Dummy 2"]],
                ]
            ],
        )
        self.assertEqual(
            famdb.resolve_names("Order"),
            [
                [
                    2,
                    True,
                    PC_HIGH,
                    [["scientific name", "Order"], ["common name", "Root Dummy 2"]],
                ],
                [
                    3,
                    False,
                    PC_HIGH,
                    [
                        ["scientific name", "Other Order"],
                        ["common name", "Root Dummy 3"],
                    ],
                ],
            ],
        )
        self.assertEqual(
            famdb.resolve_names("Other Order"),
            [
                [
                    3,
                    True,
                    PC_HIGH,
                    [
                        ["scientific name", "Other Order"],
                        ["common name", "Root Dummy 3"],
                    ],
                ]
            ],
        )

    def test_get_accessions_filtered(self):
        famdb = TestDatabase.famdb

        self.assertEqual(
            sorted(list(famdb.get_accessions_filtered())),
            [
                "DR000000001",
                "DR_Repeat1",
                "TEST0001",
                "TEST0002",
                "TEST0003",
                "TEST0004",
            ],
        )
        self.assertEqual(
            list(famdb.get_accessions_filtered(tax_id=3)),
            ["TEST0002", "TEST0003"],
        )
        self.assertEqual(
            list(famdb.get_accessions_filtered(tax_id=3, ancestors=True)),
            ["TEST0001", "TEST0002", "TEST0003"],
        )
        self.assertEqual(list(famdb.get_accessions_filtered(stage=30)), ["TEST0003"])
        self.assertEqual(list(famdb.get_accessions_filtered(stage=60)), [])
        self.assertEqual(
            list(famdb.get_accessions_filtered(is_hmm=True, stage=10)),
            [],
        )
        self.assertEqual(
            list(famdb.get_accessions_filtered(is_hmm=False, stage=10)),
            ["TEST0004"],
        )
        self.assertEqual(list(famdb.get_accessions_filtered(stage=10, is_hmm=True)), [])
        self.assertEqual(
            list(famdb.get_accessions_filtered(name="Test family TEST0004")),
            ["TEST0004"],
        )
        self.assertEqual(
            list(famdb.get_accessions_filtered(repeat_type="SINE")), ["TEST0004"]
        )

        self.assertEqual(
            list(famdb.get_accessions_filtered(tax_id=4, descendants=True)),
            ["TEST0004", "DR_Repeat1"],
        )
        # curated/uncurated are backwards because it's easier than rewriting all the family names and all the tests
        self.assertEqual(
            list(famdb.get_accessions_filtered(uncurated_only=True)),
            ["DR000000001"],
        )
        self.assertEqual(
            list(famdb.get_accessions_filtered(curated_only=True)),
            ["TEST0001", "TEST0002", "TEST0003", "TEST0004", "DR_Repeat1"],
        )
