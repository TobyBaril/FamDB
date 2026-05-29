# FamDB Quickstart

This guide covers the minimal steps to install FamDB with Dfam curated consensus
families, optionally augment with RepBase, and run common queries.

---

## Prerequisites

Python 3.6 or later with the `h5py` package:

```
pip3 install --user h5py
```

Download the release from github:

```

```


---

## Step 1: Download Dfam component files

Run the interactive downloader from the FamDB install directory:

```
python3 utils/download_dfam.py
```

The script connects to the Dfam server and presents a menu of available
components and partitions.  For most use cases, download:

- **Root** (required) -- taxonomy index, no family data
- **Curated consensus** -- curated families with consensus sequences

Select `all` partitions for each chosen component when prompted.  Files are
written to `Libraries/famdb/` by default.

> **Output directory:** Use `-o /path/to/dir` to install elsewhere, then
> pass that path to `famdb.py -i` in all subsequent commands.

---

## Step 2: Add RepBase families (optional)

RepBase RepeatMasker Edition is available from GIRI (<https://www.girinst.org/>).
After downloading, copy `RMRBSeqs.embl` into the `Libraries/` directory:

```
cp /path/to/RepBaseRepeatMaskerEdition/RMRBSeqs.embl Libraries/
```

Then merge into the installed FamDB partitions:

```
python3 utils/merge_repbase.py -i Libraries/famdb
```

The script is idempotent: re-running it is safe and skips families already
merged.

---

## Step 3: Verify the installation

```
python3 famdb.py -i Libraries/famdb info
```

Expected output:

```
[paste output here]
```

---

## Common queries

For convenience, set a shell variable for the database directory:

```
set FAMDB = Libraries/famdb          # tcsh / csh
export FAMDB=Libraries/famdb         # bash / sh
```

### Check which partitions are needed for a species

```
python3 famdb.py -i $FAMDB check 'Mus musculus'
```

```
[paste output here]
```

### Browse the taxonomy lineage with family counts

```
python3 famdb.py -i $FAMDB lineage -ad 'Mus musculus'
```

```
[paste output here]
```

### Export families in FASTA format

```
python3 famdb.py -i $FAMDB families -ad -f fasta_name 'Mus musculus'
```

```
[paste output here]
```

To restrict to a particular repeat class:

```
python3 famdb.py -i $FAMDB families -ad -f fasta_name --class LINE 'Mus musculus'
```

To export in EMBL format with full metadata:

```
python3 famdb.py -i $FAMDB families -ad -f embl 'Mus musculus'
```

---

## Further reading

- `python3 famdb.py <command> --help` -- full option list for any command
- `README.md` -- complete command reference and file format documentation
