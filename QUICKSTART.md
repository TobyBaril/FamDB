# FamDB Quickstart

This guide covers the minimal steps to install FamDB with Dfam curated consensus
families, optionally augment with RepBase, and run common queries.

---

## Prerequisites

Python 3.6 or later with the `h5py` package:

```
% pip3 install --user h5py
```

Download the release from github:
https://github.com/Dfam-consortium/FamDB/releases/latest

```
% wget https://github.com/Dfam-consortium/FamDB/archive/refs/tags/#.#.#.tar.gz
% tar zxvf #.#.#.tar.gz
```

---

## Step 1: Download Dfam component files

Run the interactive downloader from the FamDB install directory.
The script connects to the Dfam server and presents a menu of available
components and partitions.  For most use cases, download:

- **Root** (required) -- taxonomy index, no family data
- **Curated consensus** -- curated families with consensus sequences

```
% cd FamDB-#.#.#
% python3 utils/download_dfam.py
#
# download_dfam.py
#
#
# To identify the minimal set of partitions to download for a given
# species or taxon, simply download the root partition first, and then
# query the release details using:
#    ./famdb.py check <species or taxon>
#
Output directory: /u3/home/rhubley/projects/Claude/FamDB/Libraries/famdb  [default]
Release: RELEASE_DFAM_4_0
Fetching file listing from https://www.dfam.org/releases/Dfam_4.0/families/FamDB/ ...
Fetching partition size estimates...

Available components:
   1. Root (taxonomy index) [required]
        [root] -- partition 0  (58 MB)
   2. Curated consensus sequences [optional]
        [curated.consensus] -- partition 0  (27 MB)
   3. Curated profile HMMs [optional]
        [curated.hmm] -- partitions 0-1 (2 total)  (49 MB to 1.9 GB compressed per partition)
   4. Uncurated consensus sequences [optional]
        [uncurated.consensus] -- partitions 0-1 (2 total)  (258 KB to 3.8 GB compressed per partition)
   5. Uncurated profile HMMs [optional]
        [uncurated.hmm] -- partitions 0-109 (110 total)  (1.1 GB to 1.8 GB compressed per partition)

  Complete Dfam 4.0 download would be ~165.9 GB compressed.

Enter numbers to download (e.g. '1,3' or 'all'):
  (Components with multiple partitions will prompt for partition selection.)
> 1,2

  Root (taxonomy index): only 1 partition (0), selecting it.

  Curated consensus sequences: only 1 partition (0), selecting it.

Downloading 2 file(s) to: /u3/home/rhubley/projects/Claude/FamDB/Libraries/famdb

[1/2] root partition 0
    dfam40.0.h5.gz.md5: 100%  0.0 MB
    dfam40.0.h5.gz: 100%  57.9 MB
    MD5 OK: dfam40.0.h5.gz
    Decompressing dfam40.0.h5.gz ... done -> dfam40.0.h5
[2/2] curated.consensus partition 0
    dfam40.curated.consensus.0.h5.gz.md5: 100%  0.0 MB
    dfam40.curated.consensus.0.h5.gz: 100%  27.0 MB
    MD5 OK: dfam40.curated.consensus.0.h5.gz
    Decompressing dfam40.curated.consensus.0.h5.gz ... done -> dfam40.curated.consensus.0.h5

--- Summary ---
  2/2 completed successfully
```

---

## Step 2: Add RepBase families (optional)

RepBase RepeatMasker Edition (20181026) is available from GIRI (<https://www.girinst.org/>).
After downloading, copy RepBaseRepeatMaskerEdition-20181026.tar.gz to the FamDB installation
directory and run:

```
% cp RepBaseRepeatMaskerEdition-20181026.tar.gz FamDB-#.#.#
% cd FamDB-#.#.#
% tar zxvf RepBaseRepeatMaskerEdition-20181026.tar.gz
Libraries/
Libraries/RMRBSeqs.embl
Libraries/README.RMRBSeqs
```

Then merge into the installed FamDB partitions:

```
% python3 utils/merge_repbase.py
INFO: Default Libraries directory: /u3/home/rhubley/projects/Claude/FamDB/Libraries
INFO: Sourcing RMRBMeta.embl from Libraries directory: /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRBMeta.embl
INFO: Sourcing RMRBSeqs.embl from Libraries directory: /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRBSeqs.embl
INFO: Combining /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRBMeta.embl + /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRBSeqs.embl -> /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRB.embl
INFO: Reading sequences from /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRBSeqs.embl ...
INFO:   Read 49,011 sequences
INFO: Reading metadata from /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRBMeta.embl ...
INFO:   Combined 49,011 records (0 skipped - no sequence)
INFO: Sourcing RMRB_DUP.txt from Libraries directory: /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRB_DUP.txt
INFO: Merge needed: 1 new CC partition(s): dfam40.curated.consensus.0.h5
INFO: Loaded 4,119 exclusion names from /u3/home/rhubley/projects/Claude/FamDB/Libraries/RMRB_DUP.txt
INFO: Opening FamDB at /u3/home/rhubley/projects/Claude/FamDB/Libraries/famdb for writing ...
INFO:   Taxonomy lookup: 25,026 entries
INFO: Added 44,852 of 44,852 families
INFO: Rebuilding pruned taxonomy tree (499 newly-valued nodes)
INFO: Pruned Tree Updated
INFO: Finalizing files ...
INFO: Done. Added 44,852 RepBase families (RepBase version: 20181026).
```

---

## Step 3: Verify the installation

```
% python3 famdb.py info

FamDB Directory               : /u3/home/rhubley/projects/Claude/FamDB/Libraries/famdb
FamDB Creation Format Version : 3.0.0
FamDB Creation Date           : 2026-05-28 13:38:48.030055

Database : Dfam withRBRM
Version  : 4.0
Date     : 2026-05-22

Dfam - A database of transposable element (TE) sequence alignments and HMMs.

2 Partitions Present
Total consensus sequences present: 75511
Total HMMs present               : 0


Installed Components
--------------------

 Curated Consensus:
     partition 0 [dfam40.curated.consensus.0.h5]:  root  75,511 families

 Curated HMMs:
     partition 0 [dfam40.curated.hmm.0.h5]:  root       --- not present ---
     partition 1 [dfam40.curated.hmm.1.h5]:  Bilateria  --- not present ---

 Uncurated Consensus:
     partition 0 [dfam40.uncurated.consensus.0.h5]:  root       --- not present ---
     partition 1 [dfam40.uncurated.consensus.1.h5]:  Eukaryota  --- not present ---

 Uncurated HMMs:
     partition 0 [dfam40.uncurated.hmm.0.h5]:      root                           --- not present ---
     partition 1 [dfam40.uncurated.hmm.1.h5]:      Pooideae                       --- not present ---
     partition 2 [dfam40.uncurated.hmm.2.h5]:      Panicoideae                    --- not present ---
     partition 3 [dfam40.uncurated.hmm.3.h5]:      Poaceae                        --- not present ---
     partition 4 [dfam40.uncurated.hmm.4.h5]:      Petrosaviidae                  --- not present ---
...
```

---

## Common queries

### Check which partitions are needed for a species

```
% python3 famdb.py check 'Mus musculus'

Partition check for 'Mus musculus' (tax id: 10090):

  Curated Consensus    partition 0 [root]:  present
  Curated HMMs         partition 0 [root]:  MISSING  [dfam40.curated.hmm.0.h5]
                       partition 1 [Bilateria]:  MISSING  [dfam40.curated.hmm.1.h5]
  Uncurated Consensus  partition 0 [root]:  MISSING  [dfam40.uncurated.consensus.0.h5]
                       partition 1 [Eukaryota]:  MISSING  [dfam40.uncurated.consensus.1.h5]
  Uncurated HMMs       partition 0 [root]:  MISSING  [dfam40.uncurated.hmm.0.h5]
                       partition 70 [Rodentia]:  MISSING  [dfam40.uncurated.hmm.70.h5]
                       partition 73 [Eutheria]:  MISSING  [dfam40.uncurated.hmm.73.h5]
                       partition 92 [Euteleostomi]:  MISSING  [dfam40.uncurated.hmm.92.h5]
                       partition 93 [Vertebrata <vertebrates>]:  MISSING  [dfam40.uncurated.hmm.93.h5]
                       partition 94 [Deuterostomia]:  MISSING  [dfam40.uncurated.hmm.94.h5]
                       partition 95 [Bilateria]:  MISSING  [dfam40.uncurated.hmm.95.h5]
                       partition 96 [Metazoa]:  MISSING  [dfam40.uncurated.hmm.96.h5]
```

### Browse the taxonomy lineage with family counts

```
% python3 famdb.py lineage -ad 'Mus musculus'
# Format: <NCBI tax ID> <scientific name> [<# families>]
#        where counts represent curated (DF) and uncurated (DR) families

1 root [9]
└─33208 Metazoa [5]
  └─7742 Vertebrata <vertebrates> [80]
    └─117571 Euteleostomi [1]
      └─32523 Tetrapoda [19]
        └─32524 Amniota [102]
          └─40674 Mammalia [67]
            └─32525 Theria <mammals> [69]
              └─9347 Eutheria [388]
                └─1437010 Boreoeutheria [40]
                  └─314146 Euarchontoglires [44]
                    └─314147 Glires [3]
                      └─9989 Rodentia [18]
                        └─1963758 Myomorpha [17]
                          └─337687 Muroidea [59]
                            └─10066 Muridae [35]
                              └─39107 Murinae [186]
                                └─10088 Mus <genus> [221]
                                  └─10090 Mus musculus [32]
                                    ├─10091 Mus musculus castaneus [1122]
                                    ├─10092 Mus musculus domesticus [1147]
                                    ├─39442 Mus musculus musculus [1162]
                                    └─57486 Mus musculus molossinus [1261]
```

### Export families in FASTA format

```
% python3 famdb.py families --curated -ad -f fasta_name 'Mus musculus' 
>7SLRNA_short_ @Rodentia [S:35,40,50]
GCCGGGCGCGGTGGCGCGTGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGTGGGAGGAT
CGCTTGAGTCCAGGAGTTCTGGGCTGTAGTGCGCTATGCCGATCGGGTGTCCGCACTAAG
TTCGGCATCAATATGGTGACCTCCCGGGAGCGGGGGACCACCAGGTTGCCTAAGGAGGGG
TGAACCGGCCCAGGTCGGAAACGGAGCAGGTCAAAACTCCCGTGCTGATCAGTAGTGGGA
TCGCGCCTGTGAATAGCCACTGCACTCCAGCCTGAGCAACATAGCGAGACCCCGTCTCTT
AAAAAAAAAAAAAA
>B1-dID @Rodentia
AGCCGGGTGTGGTGGCGCAYGCCTGTAATCCCAGCSACTTGGGAGGCTGAGGCAGGAGGA
TCACAAGTTCAAGGCCAGCCTCAGCAACTTAGTGAGGCCCTAAGCAACTTAGTGAGACCC
TGTCTCAAAATAAAAAAAAAAAAAAAAGGGGCTGGGGATGTGGCTCAGTGGTAGAGTGCC
CCTGGGTTCAATCCCCAGTACCAAAAAAAAAAAAAAAAAAA
...
```

To restrict to a particular repeat class:

```
% python3 famdb.py families --curated -ad -f fasta_name --class LINE 'Mus musculus'
>L1M5_orf2 @Theria_mammals [S:55,75]
ATGGTAGATTTAAACCCAANCATATCAATAATTACATTAAATGTAAATGGNCTAAACACT
CCAATTAAAAGGCAGAGATTGTCAGACTGGATAAAAAAACAAGACCCAACTATATGCTGT
CTACAAGAGACGCACTTTAAATATAAAGACACAGANAGGTTGAAAGTAAAAGGATGGAAA
...
```

To export in EMBL format with full metadata:

```
% python3 famdb.py families --curated -ad -f embl 'Mus musculus'
CC   Dfam - A database of transposable element (TE) sequence alignments and HMMs
CC   Copyright (C) 2026 The Dfam consortium.
CC   
CC   Release: Dfam_4.0
CC   Date   : 2026-05-22
CC   
...
```

---

## Further reading

- `% python3 famdb.py <command> --help` -- full option list for any command
- `README.md` -- complete command reference and file format documentation
