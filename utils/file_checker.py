#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    This script can be used to check the data contained in FamDB files. 
    It performs the following tests:
        Open the file as an h5py file
        Load and display metadata including version info and expected counts
        Test for the existence of the datasets expected in FamDB files
        Attempt to use the FamDBRoot or FamDBLeaf classes to read the file
        Check for any interruptions in the file's change log.
    
    Note that this does require h5py and a FamDB installation to function.

    Usage: file_checker.py /input/path

    input    : A path to the file to check.

SEE ALSO:
    famdb.py
    Dfam: http://www.dfam.org

AUTHOR(S):
    Anthony Gray <agray@systemsbiology.org>

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

import sys
import os
import logging
import argparse
import h5py
import datetime

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from famdb_classes import FamDBLeaf, FamDBRoot
from famdb_globals import (
    FILE_DESCRIPTION,
    GROUP_FAMILIES,
    GROUP_LOOKUP_BYNAME,
    GROUP_LOOKUP_BYTAXON,
    GROUP_LOOKUP_BYSTAGE,
    GROUP_NODES,
    GROUP_REPEATPEPS,
    GROUP_FILE_HISTORY,
    DATA_PARTITION_CACHE,
    DATA_NAMES_CACHE,
    META_DB_VERSION,
    META_DB_NAME,
    META_CREATED,
    META_FAMDB_VERSION,
    FAMDB_COMPONENT_FILE_RE,
    FAMDB_ROOT_FILE_RE,
)

# Expected groups/datasets per file type in v3 FamDB.
# "root" entries must be present in root files.
# "component" entries must be present in component (leaf) files.
ROOT_REQUIRED = [
    GROUP_NODES,
    GROUP_REPEATPEPS,
    GROUP_LOOKUP_BYTAXON,
    DATA_PARTITION_CACHE,
    DATA_NAMES_CACHE,
    GROUP_FILE_HISTORY,
]

COMPONENT_REQUIRED = [
    GROUP_FAMILIES,
    GROUP_LOOKUP_BYNAME,
    GROUP_LOOKUP_BYSTAGE,
    GROUP_FILE_HISTORY,
]


def attempt(func):
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print(f"Function {func.__name__} Failed: {e}")
            exit(1)

    return wrapper


@attempt
def file_is_root(path):
    filename = os.path.basename(path)
    if FAMDB_COMPONENT_FILE_RE.match(filename):
        return False
    if FAMDB_ROOT_FILE_RE.match(filename):
        return True
    raise Exception(
        f"File '{filename}' does not match a recognized FamDB v3 filename pattern.\n"
        "Expected '<base>.0.h5' (root) or '<base>.(curated|uncurated).(consensus|hmm).<N>.h5' (component)."
    )


@attempt
def get_group(file, group):
    thing = file.get(group)
    if thing is None:
        return False
    has_keys = callable(getattr(thing, "keys", None))
    if has_keys:
        keys = list(thing.keys())
        if group == GROUP_FAMILIES:
            # All families are sorted into DF, DR, or Aux bins
            return set(keys) <= {"Aux", "DF", "DR"}
        if group == GROUP_LOOKUP_BYNAME:
            return len(keys) > 0
        if group == GROUP_LOOKUP_BYSTAGE:
            return len(keys) > 0
        if group == GROUP_LOOKUP_BYTAXON:
            return len(keys) > 0
        if group == GROUP_NODES:
            return all(x.isdigit() for x in keys)
        if group == GROUP_FILE_HISTORY:
            for elem in keys:
                try:
                    datetime.datetime.strptime(elem, "%Y-%m-%d %H:%M:%S.%f")
                except ValueError:
                    return False
            return True
    else:
        # Dataset (not a group)
        if group in (GROUP_REPEATPEPS, DATA_PARTITION_CACHE, DATA_NAMES_CACHE):
            return True
    return False


@attempt
def get_famdb_class(is_root, path):
    if is_root:
        return FamDBRoot(path, "r")
    else:
        return FamDBLeaf(path, "r")


@attempt
def get_attrs(file):
    print(
        "  Metatadata Retrieved:\n"
        f"    FamDB Version: {file.attrs[META_FAMDB_VERSION] if file.attrs.get(META_FAMDB_VERSION) else None}\n"
        f"    File Created: {file.attrs[META_CREATED]}\n"
        f"    Name: {file.attrs[META_DB_NAME]}\n"
        f"    Dfam Version: {file.attrs[META_DB_VERSION]}\n"
        f"    Consensi Count: {file.attrs['count_consensus']}\n"
        f"    HMM Count: {file.attrs['count_hmm']}"
    )


def main():
    logging.basicConfig()

    parser = argparse.ArgumentParser(
        description=FILE_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input")
    args = parser.parse_args()

    is_root_file = file_is_root(args.input)

    try:
        print(f"Testing that {args.input} can be opened as H5 file...")
        with h5py.File(args.input, "r") as file:
            get_attrs(file)
            all_expected_groups = True
            required = ROOT_REQUIRED if is_root_file else COMPONENT_REQUIRED
            for group in required:
                print(f"  Testing: {group}")
                if not get_group(file, group):
                    print(f"\t {group} Not found, or not as expected")
                    all_expected_groups = False
            # Component files must have a component_type attribute
            if not is_root_file:
                ct = file.attrs.get("component_type")
                if ct:
                    print(f"  component_type: {ct.decode() if isinstance(ct, bytes) else ct}")
                else:
                    print("\t component_type attribute missing from component file")
                    all_expected_groups = False
            # Root files must NOT have GROUP_LOOKUP_BYTAXON directly under a component key
            # (just check that GROUP_LOOKUP_BYTAXON is absent from component files)
            if not is_root_file and file.get(GROUP_LOOKUP_BYTAXON):
                print(f"\t {GROUP_LOOKUP_BYTAXON} should not be present in component files")
                all_expected_groups = False
            if all_expected_groups:
                print("  All Expected Groups Found")

    except Exception as e:
        print(f"{args.input} Could not be Opened as an H5 File: {e}")
        exit(1)

    try:
        print(
            f"Testing that {args.input} can be opened with the {'FamDBRoot' if is_root_file else 'FamDBLeaf'} class..."
        )
        with get_famdb_class(is_root_file, args.input) as file:
            if file.interrupt_check():
                print(
                    f"{args.input} Was Interrupted During Writing. Corruption Possible"
                )
            else:
                print(f"  No Interruption Detected")

    except Exception as e:
        print(f"{args.input} Could Not Be Opened As An H5 File: {e}")
        exit(1)


if __name__ == "__main__":
    main()
