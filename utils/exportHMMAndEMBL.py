#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Export the dfam database to FamDB format.

    Usage: exportHMMAndEMBL.py.py [-h] [-l LOG_LEVEL]
                         [--target_size #]
                         [--dfam_config <DFAM_CONFIG>]
                         [--cpu #][--quiet][--db-version <DB_VERSION>]
                         [--db-date <DB_DATE>]


    Export the Dfam database to HMM and EMBL format for release.
    This creates the following files:

             Dfam-#.hmm.gz
             Dfam-#.hmm.gz.md5
             Dfam-curated_only-#.hmm.gz
             Dfam-curated_only-#.hmm.gz.md5
             Dfam-#.embl.gz
             Dfam-#.embl.gz.md5
             Dfam-curated_only-#.embl.gz
             Dfam-curated_only-#.embl.gz.md5

    It does so by batching the extraction from the MySQL database up to
    cpu # of batches.  Each batch is processed in parallel.
    The final number of output files is controlled by the target_size
    parameter.  The output files are compressed using pigz and an MD5


    Args:
        --help, -h             : Show this help message and exit.
        --dfam_confg, -c       : Use a specific dfam config file rather than
                                   using the DFAM_CONF environment variable,
                                   or a relative path.
        --quiet, -q            : Reduce the amount of logging output.
        --cpu, -p              : The number of CPUs to assign to the nhmmer
                                   searches. Default: 16 cpus
        --db-version           : Set the database version explicitly, overriding
                                    values in the database.
        --db-date              : Set the database date explicitly, overriding
                                    values in the database.
        --target_size          : The target size for the output files. Default:
                                    100_000_000_000 bytes (100GB)

SEE ALSO:
    famdb.py
    Dfam: http://www.dfam.org

AUTHOR(S):
    Robert Hubley <rhubley@systemsbiology.org>
    Anthony Gray <agray@systemsbiology.org>
    Jeb Rosen <jeb.rosen@systemsbiology.org>

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
import datetime
import itertools
import json
import logging
import re
import time
import sys
import textwrap
import os
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, wait, as_completed
import traceback

# Import SQLAlchemy
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import select, func, over, cast, text, Integer
from sqlalchemy.sql import literal_column, over
from sqlalchemy.orm import aliased

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
from famdb_data_loaders import iterate_db_families

LOGGER = logging.getLogger(__name__)


##
## Imported from famdb_globals.py
##
COPYRIGHT_TEXT = """Dfam - A database of transposable element (TE) sequence alignments and HMMs
Copyright (C) %s The Dfam consortium.

Release: Dfam_%s
Date   : %s

This database is free; you can redistribute it and/or modify it
as you wish, under the terms of the CC0 1.0 license, a
'no copyright' license:

The Dfam consortium has dedicated the work to the public domain, waiving
all rights to the work worldwide under copyright law, including all related
and neighboring rights, to the extent allowed by law.

You can copy, modify, distribute and perform the work, even for
commercial purposes, all without asking permission.
See Other Information below.

Other Information

o In no way are the patent or trademark rights of any person affected by
  CC0, nor are the rights that other persons may have in the work or in how
  the work is used, such as publicity or privacy rights.
o Makes no warranties about the work, and disclaims liability for all uses of the
  work, to the fullest extent permitted by applicable law.
o When using or citing the work, you should not imply endorsement by the Dfam consortium.

You may also obtain a copy of the CC0 license here:
http://creativecommons.org/publicdomain/zero/1.0/legalcode
"""


def compress_and_signature_files(output_files, compression_level=5):
    """
    Compresses output files using pigz and generates an MD5 checksum for each file.

    Parameters:
    - output_files (list): List of output file paths to process.
    - compression_level (int): Compression level for pigz (default: 5).

    Returns:
    - List of tuples with (compressed_file, file_size, md5_file).
    """
    processed_files = []

    for file_path in output_files:
        if not os.path.exists(file_path):
            LOGGER.error(f"File not found: {file_path}")
            continue

        # Calculate original file size
        file_size = os.path.getsize(file_path)
        # LOGGER.info(f"File: {file_path} Size: {file_size} bytes")

        # Compress the file using pigz -5
        compressed_file = f"{file_path}.gz"
        try:
            # LOGGER.info(f"Compressing {file_path} to {compressed_file} using pigz -{compression_level}")
            subprocess.run(
                ["pigz", f"-{compression_level}", "-f", file_path], check=True
            )
        except subprocess.CalledProcessError as e:
            LOGGER.error(f"Compression failed for {file_path}: {e}")
            continue

        # Calculate MD5 checksum
        md5_file = f"{compressed_file}.md5"
        try:
            md5sum = calculate_md5(compressed_file)
            with open(md5_file, "w") as md5_handle:
                md5_handle.write(f"{md5sum}  {os.path.basename(compressed_file)}\n")
            # LOGGER.info(f"MD5 checksum for {compressed_file}: {md5sum}")
        except Exception as e:
            LOGGER.error(f"Failed to generate MD5 for {compressed_file}: {e}")
            continue

        processed_files.append((compressed_file, file_size, md5_file))

    return processed_files


def calculate_md5(file_path, buffer_size=1024 * 1024):
    """
    Calculate the MD5 checksum of a file.

    Parameters:
    - file_path (str): Path to the file.
    - buffer_size (int): Size of buffer for reading the file in chunks (default: 1MB).

    Returns:
    - str: The MD5 checksum in hexadecimal format.
    """
    import hashlib

    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        while True:
            chunk = f.read(buffer_size)
            if not chunk:
                break
            md5.update(chunk)
    return md5.hexdigest()


def concatenate_files(
    file_list,
    output_file_prefix,
    file_type,
    max_record_size,
    db_version,
    db_date,
    remove_originals=True,
    buffer_size=1024 * 1024,
):
    """
    Concatenate a list of files into multiple output files with a maximum size limit,
    ensuring that records are not split across files. The record terminator is "//".

    Parameters:
    - file_list (list): List of file paths to concatenate.
    - output_file_prefix (str): Prefix for output files (e.g., "output").
    - file_type (str): Type of file to generate ("hmm" or "embl").
    - max_record_size (int): Maximum size (in bytes) for each output file.
    - db_version (str): Version of the database for header generation.
    - db_date (str): Date of the database for header generation.
    - remove_originals (bool): Whether to remove original files after concatenation (default: True).
    - buffer_size (int): Size of the read buffer in bytes (default: 1MB).

    Returns:
    - List of generated output file paths.
    """

    # Validate file type and set file extension and header format
    if file_type not in ["hmm", "embl"]:
        LOGGER.error(f"Unsupported file type '{file_type}'. Use 'hmm' or 'embl'.")
        sys.exit(1)

    extension = f".{file_type}"
    header_prefix = "#   " if file_type == "hmm" else "CC   "
    header_text = re.sub(
        r"(?m)^",
        header_prefix,
        COPYRIGHT_TEXT % (datetime.datetime.now().year, db_version, db_date),
    )

    # Check if all files exist before proceeding
    missing_files = [file for file in file_list if not os.path.exists(file)]
    if missing_files:
        LOGGER.error(f"The following files do not exist:\n" + "\n".join(missing_files))
        sys.exit(1)

    output_file_index = 1
    current_size = 0
    output_handle = None
    record_buffer = []
    output_files = []

    for file_path in file_list:
        # LOGGER.info(f"Processing {file_path}")
        with open(file_path, "r") as input_handle:
            while True:
                line = input_handle.readline()
                if not line:
                    break

                record_buffer.append(line)

                # Check for the record terminator ("//")
                if line.strip() == "//":
                    record_size = sum(len(l) for l in record_buffer)

                    # Open a new output file if needed
                    if not output_handle or (
                        current_size + record_size > max_record_size
                    ):
                        if output_handle:
                            output_handle.close()

                        output_filename = (
                            f"{output_file_prefix}-{output_file_index}{extension}"
                        )
                        # LOGGER.info(f"Creating new output file: {output_filename}")
                        output_handle = open(output_filename, "w")
                        output_handle.write(header_text + "\n")
                        current_size = len(header_text) + 1  # +1 for newline
                        output_files.append(output_filename)
                        output_file_index += 1

                    # Warn if a single record exceeds max_record_size
                    # if record_size > max_record_size:
                    #    LOGGER.warning(
                    #        f"Record size ({record_size} bytes) exceeds max_record_size ({max_record_size} bytes)."
                    #        " The file size will exceed the limit."
                    #    )

                    # Write the complete record to the file
                    output_handle.writelines(record_buffer)
                    current_size += record_size
                    record_buffer = []

        # Remove the original file if required
        if remove_originals:
            # LOGGER.info(f"Removing {file_path}")
            os.remove(file_path)

    # Close the final output file handle if open
    if output_handle:
        output_handle.close()

    # LOGGER.info(f"Concatenation completed. {len(output_files)} files created with prefix '{output_file_prefix}'.")

    return output_files


##
## Imported from famdb_helper_classes and modified to get tax_db from parameters
##
def to_embl(
    family, tax_db, include_meta=True, include_seq=True
):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Converts 'family' to EMBL format."""

    if include_seq and family.consensus is None:
        # Skip families without consensus sequences, if sequences were required.
        # metadata-only formats will still include families without a consensus sequence.
        return None

    sequence = family.consensus or ""

    out = ""

    # Appends to 'out':
    # "TAG  Text"
    #
    # Or if wrap=True and 'text' has multiple lines:
    # "TAG  Line 1"
    # "TAG  Line 2"
    def append(tag, text, wrap=False):
        nonlocal out
        if not text:
            return

        prefix = "%-5s" % tag
        if wrap:
            text = textwrap.fill(str(text), width=72)
        out += textwrap.indent(str(text), prefix)
        out += "\n"

    # Appends to 'out':
    # "FT                   line 1"
    # "FT                   line 2"
    def append_featuredata(text):
        nonlocal out
        prefix = "FT                   "
        if text:
            out += textwrap.indent(textwrap.fill(str(text), width=72), prefix)
            out += "\n"

    id_line = family.accession
    if family.version is not None:
        id_line += "; SV " + str(family.version)

    append("ID", "%s; linear; DNA; STD; UNC; %d BP." % (id_line, len(sequence)))
    append("NM", family.name)
    out += "XX\n"
    append("AC", family.accession + ";")
    out += "XX\n"
    append("DE", family.title, True)
    out += "XX\n"

    if include_meta:
        if family.aliases:
            for alias_line in family.aliases.splitlines():
                [db_id, db_link] = map(str.strip, alias_line.split(":"))
                if db_id == "Repbase":
                    append("DR", "Repbase; %s." % db_link)
                    out += "XX\n"

        if family.repeat_type == "LTR":
            append(
                "KW",
                "Long terminal repeat of retrovirus-like element; %s." % family.name,
            )
        else:
            append(
                "KW", "%s/%s." % (family.repeat_type or "", family.repeat_subtype or "")
            )
        out += "XX\n"

        for clade_id in family.clades:
            lineage_str = tax_db[clade_id][1]
            lineage = lineage_str.split(";")
            if lineage[0] == "root":
                lineage = lineage[1:]
            if len(lineage) > 0:
                append("OS", lineage[-1])
                append("OC", "; ".join(lineage[:-1]) + ".", True)
        out += "XX\n"

        if family.citations:
            citations = json.loads(family.citations)
            citations.sort(key=lambda c: c["order_added"])
            for cit in citations:
                append(
                    "RN", "[%d] (bases 1 to %d)" % (cit["order_added"], family.length)
                )
                append("RA", cit["authors"], True)
                append("RT", cit["title"], True)
                append("RL", cit["journal"])
                out += "XX\n"

        append("CC", family.description, True)
        out += "CC\n"
        append("CC", "RepeatMasker Annotations:")
        append("CC", "     Type: %s" % (family.repeat_type or ""))
        append("CC", "     SubType: %s" % (family.repeat_subtype or ""))

        species_names = [tax_db[c][0] for c in family.clades]
        append("CC", "     Species: %s" % ", ".join(species_names))

        append("CC", "     SearchStages: %s" % (family.search_stages or ""))
        append("CC", "     BufferStages: %s" % (family.buffer_stages or ""))
        if family.refineable:
            append("CC", "     Refineable")

        if family.coding_sequences:
            out += "XX\n"
            append("FH", "Key             Location/Qualifiers")
            out += "FH\n"
            for cds in json.loads(family.coding_sequences):
                for element in cds:
                    if type(cds[element]) == str:
                        cds[element] = cds[element].replace('"', "")

                append(
                    "FT",
                    "CDS             %d..%d" % (cds["cds_start"], cds["cds_end"]),
                )
                append_featuredata('/product="%s"' % cds["product"])
                append_featuredata("/number=%s" % cds["exon_count"])
                append_featuredata('/note="%s"' % cds["description"])
                append_featuredata('/translation="%s"' % cds["translation"])

        out += "XX\n"

    if include_seq:
        sequence = sequence.lower()
        i = 0
        counts = {"a": 0, "c": 0, "g": 0, "t": 0, "other": 0}
        for char in sequence:
            if char not in counts:
                char = "other"
            counts[char] += 1

        append(
            "SQ",
            "Sequence %d BP; %d A; %d C; %d G; %d T; %d other;"
            % (
                len(sequence),
                counts["a"],
                counts["c"],
                counts["g"],
                counts["t"],
                counts["other"],
            ),
        )

        while i < len(sequence):
            chunk = sequence[i : i + 60]
            i += 60

            j = 0
            line = ""
            while j < len(chunk):
                line += chunk[j : j + 10] + " "
                j += 10

            out += "     %-66s %d\n" % (line, min(i, len(sequence)))

    out += "//\n"

    return out


def to_dfam_hmm(family, tax_db):  # pylint: disable=too-many-locals,too-many-branches
    """
    Converts 'family' to Dfam-style HMM format.
    'famdb' is used for lookups in the taxonomy database (id -> name).
    """
    if family.model is None:
        return None

    out = ""

    # Appends to 'out':
    # "TAG   Text"
    #
    # Or if wrap=True and 'text' has multiple lines:
    # "TAG   Line 1"
    # "TAG   Line 2"
    def append(tag, text, wrap=False):
        nonlocal out
        if not text:
            return

        prefix = "%-6s" % tag
        text = str(text)
        if wrap:
            text = textwrap.fill(text, width=72)
        out += textwrap.indent(text, prefix)
        out += "\n"

    # TODO: Compare to e.g. finditer(). This does a lot of unnecessary
    # allocation since most of model_lines are appended verbatim.
    model_lines = family.model.split("\n")

    i = 0
    for i, line in enumerate(model_lines):
        if line.startswith("HMMER3"):
            out += line + "\n"

            name = family.name or family.accession
            append("NAME", name)
            append("ACC", family.accession_with_optional_version())
            append("DESC", family.title)
        elif any(map(line.startswith, ["NAME", "ACC", "DESC"])):
            # Correct version of this line was output already
            pass
        elif line.startswith("CKSUM"):
            out += line + "\n"
            break
        else:
            out += line + "\n"

    th_lines = []
    if family.taxa_thresholds:
        for threshold in family.taxa_thresholds.split("\n"):
            parts = threshold.split(",")
            tax_id = int(parts[0])
            try:
                (hmm_ga, hmm_tc, hmm_nc, hmm_fdr) = map(float, parts[1:])
            except Exception as err:
                hmm_ga = 0.0
                hmm_tc = 0.0
                hmm_nc = 0.0
                hmm_fdr = 0.0
                print(
                    f"Error in thresholds for accession={family.accession_with_optional_version()} and taxid={tax_id}",
                    file=sys.stderr,
                )

            # only recover name, do need for partition number
            tax_name = tax_db[tax_id][0]
            th_lines += [
                "TaxId:%d; TaxName:%s; GA:%.2f; TC:%.2f; NC:%.2f; fdr:%.3f;"
                % (tax_id, tax_name, hmm_ga, hmm_tc, hmm_nc, hmm_fdr)
            ]

    if family.general_cutoff:
        append("GA", "%.2f;" % family.general_cutoff)
        append("TC", "%.2f;" % family.general_cutoff)
        append("NC", "%.2f;" % family.general_cutoff)

    for th_line in th_lines:
        append("TH", th_line)

    if family.build_method:
        append("BM", family.build_method)
    if family.search_method:
        append("SM", family.search_method)

    append("CT", (family.classification and family.classification.replace("root;", "")))

    for clade_id in family.clades:
        tax_name = tax_db[clade_id][0]
        append("MS", "TaxId:%d TaxName:%s" % (clade_id, tax_name))

    append("CC", family.description, True)
    append("CC", "RepeatMasker Annotations:")
    append("CC", "     Type: %s" % (family.repeat_type or ""))
    append("CC", "     SubType: %s" % (family.repeat_subtype or ""))

    species_names = [tax_db[c][0] for c in family.clades]
    append("CC", "     Species: %s" % ", ".join(species_names))

    append("CC", "     SearchStages: %s" % (family.search_stages or ""))
    append("CC", "     BufferStages: %s" % (family.buffer_stages or ""))

    if family.refineable:
        append("CC", "     Refineable")

    # Append all remaining lines unchanged
    out += "\n".join(model_lines[i + 1 :])

    return out


def export_families(
    args, session, tax_db, start, end
):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    """Exports from a Dfam database to a FamDB file."""

    to_import = []
    target_count = 0

    limit = None
    query = (
        session.query(dfam.Family).filter(dfam.Family.id.between(start, end))
    ).limit(limit)

    # if not args.include_uncurated and not args.db_partition:
    #    query = query.filter(dfam.Family.accession.like("DF%"))
    # query = query.filter(dfam.Family.disabled != 1)
    hmm_curated_file = open(
        "batch-" + str(start) + "-" + str(end) + "-curated_only.hmm", "w"
    )
    embl_curated_file = open(
        "batch-" + str(start) + "-" + str(end) + "-curated_only.embl", "w"
    )
    hmm_file = open("batch-" + str(start) + "-" + str(end) + ".hmm", "w")
    embl_file = open("batch-" + str(start) + "-" + str(end) + ".embl", "w")

    target_count += query.count()
    # LOGGER.info(f"Including {target_count} families from database")

    to_import = itertools.chain(to_import, iterate_db_families(session, query))

    start = time.perf_counter()
    report_start = start
    # Note about timing.  At this stage we haven't executed the iterate_db_families function yet
    # to iterate over the yielded family objects.  Therefore, the first time through this loop there
    # will be some overhead while it loads the classification nodes.  The remaining cycles will only
    # include the inner yeild loop in iterate_db_families.
    report_every = 1000
    if target_count > 1000000:
        report_every = int(target_count / 10000)

    count = 0
    for family in to_import:
        count += 1

        hmm_data = to_dfam_hmm(family, tax_db)
        # print("HMM RECORD:\n" + foo)
        embl_data = to_embl(family, tax_db)
        # print("EMBL RECORD:\n" + foo)

        hmm_file.write(hmm_data)
        embl_file.write(embl_data)
        if family.accession.startswith("DF"):
            hmm_curated_file.write(hmm_data)
            embl_curated_file.write(embl_data)

        # LOGGER.debug(
        #    f"Added family {family.name} ({family.accession})"
        # )

        if (count % report_every) == 0:
            current = time.perf_counter()
            total_elapsed = current - start
            report_elapsed = current - report_start
            avg_time_per = total_elapsed / count
            curr_time_per = report_elapsed / report_every
            rem_time = str(
                datetime.timedelta(seconds=(curr_time_per * (target_count - count)))
            )
            LOGGER.info(
                f"Stat: {count:5d} / {target_count:5d} : {avg_time_per:.3f} avg secs per family :{curr_time_per:.3f} curr secs per family : {rem_time} HH:MM:SS remaining"
            )
            report_start = time.perf_counter()

    delta = str(datetime.timedelta(seconds=time.perf_counter() - start))
    LOGGER.info(f"Added {count} families in {delta}")


def run_export(
    args,
    conf,
    start,
    end,
    tax_db,
):
    LOGGER.info(f"\tExporting chunk {start} - {end}")

    df_engine = create_engine(conf.getDBConnStrWPassFallback("Dfam"))
    df_sfactory = sessionmaker(df_engine)
    session = df_sfactory()

    export_families(args, session, tax_db, start, end)


def main():
    """Parses command-line arguments and runs the import."""

    logging.basicConfig()

    parser = argparse.ArgumentParser(
        description=(
            "Export Dfam family data as HMM and EMBL files, split into partitions "
            "for use in FamDB releases. Connects to a Dfam database using the Dfam "
            "config file and processes families in parallel batches."
        )
    )
    parser.add_argument(
        "-l", "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR). Default: %(default)s",
    )
    parser.add_argument(
        "-t", "--target_size",
        dest="target_size",
        type=int,
        default=100_000_000_000,
        help=(
            "Target size in bytes for each output partition file. "
            "Default: %(default)s (100 GB)"
        ),
    )
    parser.add_argument(
        "-c", "--dfam_config",
        dest="dfam_config",
        help=(
            "Path to the Dfam config file. If not specified, falls back to the "
            "DFAM_CONF environment variable, then ../Conf/dfam.conf."
        ),
    )
    parser.add_argument(
        "-p", "--cpu",
        type=int,
        default=16,
        dest="cpu",
        help="Number of parallel worker processes to use. Default: %(default)s",
    )
    parser.add_argument(
        "-q", "--quiet",
        dest="quiet",
        action="store_true",
        help="Suppress the parameter summary printed at startup. Default: off",
    )
    parser.add_argument(
        "--db-version",
        help=(
            "Override the Dfam database version string written into output file "
            "headers. If not specified, the value is read from the database."
        ),
    )
    parser.add_argument(
        "--db-date",
        help=(
            "Override the Dfam database release date written into output file "
            "headers (format: YYYY-MM-DD). If not specified, the value is read "
            "from the database."
        ),
    )

    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    #   Search order: --dfam_config <path>, environment DFAM_CONF,
    #                 and finally "../Conf/dfam.conf"
    #
    conf = dc.DfamConfig(args.dfam_config)

    df_engine = create_engine(conf.getDBConnStrWPassFallback("Dfam"))
    df_sfactory = sessionmaker(df_engine)
    session = df_sfactory()
    df_ver = dfVersion.DfamVersion()
    version = df_ver.version_string

    version_info = session.query(dfam.DbVersion).one()
    db_version = version_info.dfam_version
    db_date = version_info.dfam_release_date.strftime("%Y-%m-%d")

    if args.db_version:
        db_version = args.db_version
    if args.db_date:
        db_date = args.db_date

    if not args.quiet:
        print(f"#\n# exportHMMAndEMBL.py {version}\n#")
        print("# !!!! VERIFY THESE VALUES: They will appear the top of each file !!!!")
        print(f"# db_version   : {db_version}")
        print(f"# db_date      : {db_date}")
        print(f"# target size  : {args.target_size}")

    total_families = session.query(func.count(dfam.Family.id)).scalar()
    threads = args.cpu
    batch_size = total_families // threads

    if not args.quiet:
        print("f# cpus            : {threads}")
        print(f"# total families  : {total_families}")
        print(f"# proccess size   : {batch_size}\n")

    # Obtain family.id starting and ending values for each batch
    session.execute(text("SET @row_number = 0;"))
    raw_query = text(
        """
    SELECT
        MIN(id) AS batch_start,
        MAX(id) AS batch_end
    FROM (
        SELECT
            id,
            CEIL((@row_number := @row_number + 1) / :batch_size) AS batch_number
        FROM family
        ORDER BY id
    ) subquery
    GROUP BY batch_number;
    """
    )

    batches = session.execute(raw_query, {"batch_size": batch_size}).fetchall()
    family_clades = session.query(
        dfam.DfamTaxdb.tax_id, dfam.DfamTaxdb.sanitized_name, dfam.DfamTaxdb.lineage
    ).all()
    tax_db = {row.tax_id: [row.sanitized_name, row.lineage] for row in family_clades}

    # for row in batches:
    #    print("Batch {} - {}".format(row.batch_start, row.batch_end))

    all_complete = True
    with ProcessPoolExecutor() as executor:
        futures = []
        partition_map = {}
        for row in batches:
            future = executor.submit(
                run_export,
                args,
                conf,
                row.batch_start,
                row.batch_end,
                tax_db,
            )
            partition_map[future] = row.batch_start
            futures.append(future)

        executor.shutdown(wait=True)
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                all_complete = False
                LOGGER.info(
                    f"Error in Batch {partition_map[future]} - {traceback.format_exc()}"
                )

    if all_complete:
        LOGGER.info("Finished export, proceeding with post-processing...")
    else:
        LOGGER.error("Errors Encountered Creating One Or More Files.")

    for suffix in ["embl", "hmm"]:
        for section in ["", "-curated_only"]:
            batch_files = []
            for row in batches:
                batch_file = (
                    "batch-"
                    + str(row.batch_start)
                    + "-"
                    + str(row.batch_end)
                    + section
                    + "."
                    + suffix
                )
                batch_files.append(batch_file)
            output_file = f"Dfam{section}"
            out_files = concatenate_files(
                batch_files,
                output_file,
                suffix,
                args.target_size,
                db_version,
                db_date,
            )
            processed_files = compress_and_signature_files(out_files)
            for compressed_file, file_size, md5_file in processed_files:
                LOGGER.info(
                    f"Processed File: {compressed_file}, Original Size: {file_size}"
                )


if __name__ == "__main__":
    main()
