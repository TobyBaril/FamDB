# FamDB

## Overview

FamDB is a modular HDF5-based export format and query tool developed for offline access
to the [Dfam] database of transposable element and repetitive DNA families.
FamDB stores family sequence models (profile HMMs and consensus sequences),
along with metadata including:

* Family names, aliases, and description
* Classification
* Taxa
* Citations and attribution

In addition, FamDB stores a subset of the NCBI Taxonomy relevant to the family
taxa represented in the files, facilitating quick extraction of
species/clade-specific family libraries. The query tool provides options for
exporting search results in a variety of common formats including EMBL, FASTA,
and HMMER HMM format. FamDB is intended for use as a read-only data store by
tools such as [RepeatMasker] as an alternative to unindexed EMBL or HMM files.

[Dfam]: https://www.dfam.org/
[RepeatMasker]: http://www.repeatmasker.org/

## File Format (v3)

Version 3 organizes families into **components** by curation status and model
type. Each component is independently partitioned across the taxonomy tree,
allowing users to install only the data relevant to their use case.

The four component types are:

| Code | Description |
|:----:|:------------|
| `cc` | Curated Consensus -- curated families with consensus sequences |
| `ch` | Curated HMMs -- curated families with profile HMMs |
| `uc` | Uncurated Consensus -- uncurated (DR-accession) families with consensus sequences |
| `uh` | Uncurated HMMs -- uncurated families with profile HMMs |

A complete installation consists of a single root file plus one file per
component partition:

```
<base>.0.h5                          root (taxonomy + index, no family data)
<base>.curated.consensus.0.h5        cc component, partition 0
<base>.curated.hmm.1.h5             ch component, partition 1
<base>.curated.hmm.2.h5             ch component, partition 2
<base>.uncurated.consensus.1.h5     uc component, partition 1
<base>.uncurated.hmm.1.h5           uh component, partition 1
...
```

The root file is always required. Component files are optional -- install only
the components needed for your use case. For example, a tool that uses only
consensus sequences needs only the `cc` and `uc` files.

All files from the same export must reside in the same directory. FamDB reports
a warning if files from different export runs are detected. Pass the directory
path to `famdb.py` via the `-i` option.

The `info` subcommand shows which components and partitions are installed and
which are available but not yet downloaded. The `check` subcommand reports
which specific partition files are required for a given species query and
whether each is locally present.

## Installation/Setup

### Dependencies

* Python 3.6 or later
* [`h5py`] for reading and writing HDF5 files

    ```
    pip3 install --user h5py
    ```

[`h5py`]: https://pypi.org/project/h5py/

### famdb.py

RepeatMasker includes a compatible version of `famdb.py`. This file should
generally not be installed or upgraded manually.

FamDB can also be downloaded separately. The latest release is at:
<https://github.com/Dfam-consortium/FamDB/releases/latest>

### Obtaining FamDB files

FamDB files for the current Dfam release are available at:
<https://www.dfam.org/releases/current/families/FamDB/>

Download the root file and the component partition files for the components
you need, placing all files in the same directory. For most RepeatMasker use
cases the curated consensus (`cc`) files are sufficient. Add curated HMM (`ch`)
files for higher sensitivity searches. Uncurated components (`uc`, `uh`) provide
the broader DR-accession content from Dfam.

## Usage

```
famdb.py -i <directory> <command> [options]
```

For full option details on any command:

```
famdb.py <command> --help
```

### Global options

| Option | Description |
|:-------|:------------|
| `-i DB_DIR` | Directory containing FamDB files (required) |
| `-e <component>` | Exclude a component type (`cc`, `ch`, `uc`, `uh`, or a comma-separated list) |
| `-l LOG_LEVEL` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |

The `-e` option is useful for testing or for temporarily working with a subset
of installed components without removing files.

### Taxonomy search

Most commands (`names`, `lineage`, `check`, `families`) accept a taxonomy
`term` argument. The term may be:

* An NCBI taxonomy identifier (e.g. `9606`)
* A full or partial scientific name (e.g. `'Homo sapiens'` or `homo`)
* A common name (e.g. `human`)

Multiple words can be given as separate arguments and are joined as a single
search string (`famdb.py names homo sapiens` is equivalent to
`famdb.py names 'homo sapiens'`).

Searches distinguish **exact matches** from **non-exact matches**. Commands
that operate on a single taxon (`lineage`, `families`, `check`) require
exactly one unambiguous result -- either a single exact match or, when there
are no exact matches, a single partial match. If the term is ambiguous or
not found, a list of candidates (or similarly-sounding alternatives) is shown
instead.

### info

Display metadata about the installed FamDB files, including the Dfam release
version, family counts per component, and which partitions are installed or
missing.

```
famdb.py -i DB_DIR info [--history]
```

The `--history` flag appends a timestamped changelog for each installed file,
showing every operation applied since creation (exports, appends, metadata
patches).

Example output:

```
FamDB Directory               : /data/famdb
FamDB Creation Format Version : 3.0.0
FamDB Creation Date           : 2025-01-15

Database : Dfam
Version  : 3.9
Date     : 2025-01-10

Dfam - A database of transposable element (TE) sequence alignments and HMMs.

Installed Components
--------------------

 Curated Consensus:
     partition 0 [dfam.curated.consensus.0.h5]:  root          125,432 families

 Curated HMMs:
     partition 0 [dfam.curated.hmm.0.h5]:        root           98,210 families
     partition 1 [dfam.curated.hmm.1.h5]:        Bilateria      31,047 families

 Uncurated Consensus:
     partition 0 [dfam.uncurated.consensus.0.h5]: root         542,100 families

 Uncurated HMMs:
     [ Not Installed ]
```

Partitions that are listed in the file map but not present on disk are shown
with `--- not present ---` in place of a family count, indicating they can be
downloaded and added to the directory.

### names

Search for taxonomy nodes by name or NCBI taxonomy ID.

```
famdb.py -i DB_DIR names [--format pretty|json] <term> [<term> ...]
```

Exact matches are listed before non-exact matches. Each result shows all known
names for the taxon (scientific name, common names, synonyms, etc.) along with
the partition key indicating which component file holds its families.

If no match is found, similarly-sounding names are suggested using soundex
matching.

The `json` format is intended for parsing by scripts; the `pretty` format is
human-readable but not reliably parseable.

Example:

```
$ famdb.py -i ./dfam names rattus

Exact Matches
=============
Taxon: 10114, Partition: cc:0,ch:0, Names: Rattus (scientific name), ...

Non-exact Matches
=================
Taxon: 10116, Partition: cc:0,ch:0, Names: Rattus norvegicus (scientific name), ...
Taxon: 10117, Partition: cc:0,ch:0, Names: Rattus rattus (scientific name), ...
```

### lineage

Display the taxonomy tree for a clade, with the number of families assigned
to each node.

```
famdb.py -i DB_DIR lineage [-a] [-d] [-k] [-c] [-u] [-f pretty|semicolon|totals] <term>
```

| Option | Description |
|:-------|:------------|
| `-a`, `--ancestors` | Include ancestor nodes up to the root |
| `-d`, `--descendants` | Include all descendant nodes |
| `-k`, `--complete` | Include nodes that have no assigned families (skipped by default) |
| `-c`, `--curated` | Count only curated families (DF accessions) |
| `-u`, `--uncurated` | Count only uncurated families (DR accessions) |
| `-f`, `--format` | Output format: `pretty` (default), `semicolon`, or `totals` |

By default the tree skips nodes with no family data. Use `-k`/`--complete`
to show every intermediate node in the full NCBI taxonomy, even those with
no directly assigned families.

The `pretty` format includes a header explaining the component partition codes
and notes that family counts reflect the full Dfam release -- locally missing
partitions are not subtracted. Use `famdb.py check` to verify local
installation status for a given species.

The `semicolon` format always implies `--ancestors` and `--complete`, producing
a full colon-delimited lineage path per matched taxon. This is suitable for
script consumption and is the format used by RepeatMasker internally.

The `totals` format prints a single summary line showing the number of families
found in ancestors versus lineage-specific entries for the queried taxon.

Examples:

```
famdb.py -i DB_DIR lineage -ad 'Homo sapiens'
famdb.py -i DB_DIR lineage -ad --format totals 9606
famdb.py -i DB_DIR lineage -f semicolon rattus
famdb.py -i DB_DIR lineage -adk 'Mus musculus'
```

### check

Report which component partition files are needed for a given species query
and whether each is locally installed.

```
famdb.py -i DB_DIR check [--component <cc|ch|uc|uh>] <term>
```

The check covers the queried taxon and all its ancestors, since ancestor-level
partitions contribute families to any query for a descendant species. For
example, a search against *Homo sapiens* needs families assigned at the
Eutheria level, the Vertebrata level, and so on, in addition to those assigned
at the species level itself.

The `--component` option may be repeated to restrict the check to specific
component types (e.g. `--component cc --component ch`).

Example:

```
$ famdb.py -i DB_DIR check 'Homo sapiens'

Partition check for 'Homo sapiens' (tax id: 9606):

  Curated Consensus    partition 0 [root]:             present
  Curated HMMs         partition 0 [root]:             present
                       partition 1 [Bilateria]:        present
  Uncurated Consensus  partition 0 [root]:             present
  Uncurated HMMs       partition 0 [root]:             present
                       partition 50 [Eutheria]:        MISSING  [dfam.uncurated.hmm.50.h5]
```

### families

Export all families for a clade, with optional filters.

```
famdb.py -i DB_DIR families [-a] [-d] [-c] [-u] [-f <format>]
    [--stage N] [--class TYPE] [--name PREFIX]
    [--add-reverse-complement] [--include-class-in-name]
    [--require-general-threshold]
    <term>
```

Without `-a` or `-d`, only families directly assigned to the named clade are
returned. Combining `-a` and `-d` returns the full set of families applicable
to that clade: those from ancestor nodes (shared with related species) plus
those specific to any descendant.

| Option | Description |
|:-------|:------------|
| `-a`, `--ancestors` | Include families from ancestor nodes |
| `-d`, `--descendants` | Include families from descendant nodes |
| `-c`, `--curated` | Return only curated families (DF accessions) |
| `-u`, `--uncurated` | Return only uncurated families (DR accessions) |
| `-f`, `--format` | Output format (see below) |
| `--stage N` | Include only families searched at RepeatMasker stage N (use `0` for families with no stage defined) |
| `--class TYPE` | Include only families with the given repeat type or type/subtype (e.g. `LTR` or `DNA/CMC`) |
| `--name PREFIX` | Include only families whose name starts with PREFIX |
| `--add-reverse-complement` | Append a reverse-complemented copy of each family (fasta formats only; used by RepeatMasker) |
| `--include-class-in-name` | Append the RepeatMasker type/subtype to the family name, e.g. `HERV16#LTR/ERVL` (hmm and fasta formats) |
| `--require-general-threshold` | Skip families that lack general score thresholds |

Search and buffer stages are a RepeatMasker concept. Each family is associated
with one or more search stages (the rounds of masking in which it is applied)
and optional buffer stages (additional rounds where it contributes to overlap
buffering). Stage 0 matches families with no stage annotation.

Supported formats:

| Format | Description |
|:-------|:------------|
| `summary` | (default) Human-readable: accession, name, classification, length |
| `hmm` | HMMER HMM profile with RepeatMasker metadata |
| `hmm_species` | Same as `hmm`, with species-specific GA/TC/NC thresholds substituted |
| `fasta_name` | FASTA with header `>MIR @Mammalia [S:40,60,65]` |
| `fasta_acc` | FASTA with header `>DF0000001.4 @Mammalia [S:40,60,65]` |
| `embl` | EMBL with full metadata and consensus sequence |
| `embl_meta` | EMBL with metadata only (no sequence) |
| `embl_seq` | EMBL with sequence only (no metadata) |

Examples:

```
famdb.py -i DB_DIR families -f embl_meta -ad --curated 'Drosophila melanogaster'
famdb.py -i DB_DIR families -f hmm -ad --curated --class LTR 7227
famdb.py -i DB_DIR families -f fasta_acc --name SVA --include-class-in-name hominid
famdb.py -i DB_DIR families --stage 40 -ad 'Mus musculus'
```

### family

Export a single family by accession or name.

```
famdb.py -i DB_DIR family [-f <format>] <accession>
```

The accession may be a Dfam accession number (e.g. `DF000000001`) or a family
name (e.g. `MIR3`). Supported formats are the same as for `families`, except
`hmm_species` is not available since no species context is provided.

Examples:

```
famdb.py -i DB_DIR family MIR3
famdb.py -i DB_DIR family --format fasta_acc DF000000001
famdb.py -i DB_DIR family --format embl MIR3
```

## HDF5 File Structure

This section describes the internal layout of FamDB v3 HDF5 files for
developers and advanced users.

### Overview

FamDB v3 uses a multi-file layout. All files in a set share the same
`uuid`, `db_version`, and `db_date` stored in both HDF5 attributes and a
`file_info` JSON blob; mismatched values cause a startup error.

File locking is disabled for read-only opens since it is unreliable on
network filesystems and unnecessary in the absence of concurrent writers.

### Root file (`<base>.0.h5`)

The root file is the entry point for all queries. It contains the full
taxonomy and lookup indexes but no family sequence data.

**HDF5 file-level attributes:**

| Attribute | Type | Description |
|:----------|:-----|:------------|
| `famdb_version` | str | Format version string, e.g. `"3.0.0"` |
| `created` | str | ISO timestamp of file creation |
| `db_name` | str | Database name, e.g. `"Dfam"` |
| `db_version` | str | Dfam release version, e.g. `"3.9"` |
| `db_date` | str | Release date (YYYY-MM-DD) |
| `db_copyright` | str | Copyright notice |
| `db_description` | str | Release description |
| `file_info` | str | JSON blob (see below) |
| `partition_num` | str | `"0"` for the root file |
| `root` | bool | `True` |
| `count_consensus` | int | Number of consensus sequences in this file |
| `count_hmm` | int | Number of HMM profiles in this file |

**`file_info` JSON schema:**

```json
{
  "meta": {
    "uuid":       "<shared UUID for this export set>",
    "db_version": "<Dfam version>",
    "db_date":    "<YYYY-MM-DD>"
  },
  "file_map": {
    "0":    { "filename": "...", "T_root": 1, "T_root_name": "root", "F_roots": [], "F_roots_names": [] },
    "cc.0": { "filename": "...", "T_root": 1, "T_root_name": "root", "F_roots": [1], "F_roots_names": [] },
    "ch.1": { "filename": "...", "T_root": 1, "T_root_name": "root", "F_roots": [1], "F_roots_names": [] },
    "uc.1": { "filename": "...", "T_root": 5, "T_root_name": "Bilateria", "F_roots": [5], "F_roots_names": [] }
  }
}
```

`file_map` keys are `"0"` for the root file and `"<component>.<N>"` for each
component partition. `T_root` is the NCBI taxon ID of the highest-level taxon
whose families are stored in that partition file.

**Groups and datasets:**

```
Taxonomy/
  <tax_id>/
    Children       int64[]   all child taxon IDs (full NCBI tree)
    Parent         int64[1]  parent taxon ID (full NCBI tree)
    Val_Children   int64[]   child IDs that have associated family data
    Val_Parent     int64[1]  nearest ancestor with family data
    TaxaNames      str[][]   [[name_class, name_value], ...] pairs

RepeatPeps         str[1]    FASTA protein sequences (for RepeatModeler)

Lookup/
  ByTaxon/
    <tax_id>/
      accessions   str[]     family accessions assigned to this taxon

PartitionCache     str[1]    JSON: {tax_id: {cc: N|null, ch: N|null, uc: N|null, uh: N|null}}
NamesCache         str[1]    JSON: {tax_id: [[name_class, name_value], ...]}

File_History/
  <YYYY-MM-DD HH:MM:SS.f>/  attributes: operation description
```

`Val_Children` / `Val_Parent` form a sparse tree that skips taxonomy nodes
with no associated family data. Most lineage traversals use this pruned tree
for performance; the `--complete` flag switches to the full `Children` /
`Parent` tree.

`PartitionCache` is loaded entirely into memory at startup to enable fast
taxon-to-partition routing without per-node HDF5 reads.

### Component files (`<base>.<curated|uncurated>.<consensus|hmm>.<N>.h5`)

Component files store the actual family data for one component type and one
partition of the taxonomy tree.

**HDF5 file-level attributes:**

Same as the root file, plus:

| Attribute | Type | Description |
|:----------|:-----|:------------|
| `component_type` | str | One of `cc`, `ch`, `uc`, `uh` |
| `partition_num` | str | Component key, e.g. `"ch.1"` |
| `root` | bool | `False` |

**Groups and datasets:**

```
Families/
  DF/                         curated families (DF accessions), binned by prefix
    <XX>/
      <accession>             dataset (0-length placeholder); metadata as HDF5 attrs
      <accession>.model       uint8[] gzip-compressed HMM bytes (HMM files only)
  DR/                         uncurated families (DR accessions), binned by prefix
    ...
  Aux/                        auxiliary families
    ...

Lookup/
  ByName/
    <family_name>             SoftLink -> /Families/<bin>/<accession>
  ByStage/
    <stage>/
      <accession>             SoftLink -> /Families/<bin>/<accession>

File_History/
  <YYYY-MM-DD HH:MM:SS.f>/   attributes: operation description
```

Families are binned into two-character prefix groups within `DF/` or `DR/`
to avoid the HDF5 performance degradation that occurs when a single group
exceeds ~500k entries.

### Family dataset attributes

Each family is stored as a zero-length HDF5 dataset with all metadata
as dataset-level attributes. Fields not applicable to the component type are
omitted (e.g. consensus fields are absent from HMM files and vice versa).

**Fields present in all component files:**

| Field | Type | Description |
|:------|:-----|:------------|
| `name` | str | Family name (e.g. `MIR`) |
| `accession` | str | Dfam accession (e.g. `DF000000001`) |
| `version` | int | Family version number |
| `length` | int | Consensus or model length in bp |
| `classification` | str | Semicolon-delimited classification path |
| `repeat_type` | str | RepeatMasker type (e.g. `SINE`) |
| `repeat_subtype` | str | RepeatMasker subtype (e.g. `MIR`) |
| `clades` | list | NCBI taxon IDs this family is assigned to |
| `search_stages` | str | Comma-separated RepeatMasker search stage numbers |
| `buffer_stages` | str | Comma-separated RepeatMasker buffer stage numbers |
| `aliases` | str | Alternative names / cross-references |
| `citations` | str | Literature references |
| `description` | str | Free-text description |
| `date_created` | str | Creation date |
| `date_modified` | str | Last modification date |

**Consensus-only fields (present in `cc` and `uc` files):**

| Field | Type | Description |
|:------|:-----|:------------|
| `consensus` | str | Consensus nucleotide sequence |

**HMM-only fields (present in `ch` and `uh` files):**

| Field | Type | Description |
|:------|:-----|:------------|
| `<acc>.model` | uint8[] | Gzip-compressed HMMER profile (sibling dataset) |
| `max_length` | int | Maximum target length for the HMM |
| `is_model_masked` | bool | Whether the model has been masked |
| `seed_count` | int | Number of seed sequences used to build the HMM |
| `build_method` | str | Tool and parameters used to build the HMM |
| `search_method` | str | Recommended search parameters |
| `taxa_thresholds` | str | Species-specific GA/TC/NC score thresholds (TH lines) |
| `general_cutoff` | float | General gathering threshold score |

The HMM model is stored as a gzip-compressed `uint8` dataset named
`<accession>.model` as a sibling to the family dataset, rather than as an
attribute. This means the model bytes are never decompressed during metadata
queries -- decompression only occurs when the model content is explicitly
requested.
