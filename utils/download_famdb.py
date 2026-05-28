#!/usr/bin/env python3
"""
download_famdb.py - Interactive downloader for FamDB component files from Dfam.

Fetches the available release from the Dfam server, lets you choose which
components and partitions to download, validates MD5 checksums, and
decompresses the files into Libraries/famdb/ (or a user-specified directory).

Recoverable: already-decompressed files are skipped; partially-downloaded
.gz files are re-verified before re-downloading.

Usage: download_famdb.py [-h] [-o OUTPUT_DIR] [-u URL] [--dry-run]
"""

import argparse
import gzip
import hashlib
import os
import re
import shutil
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from html.parser import HTMLParser

DEFAULT_URL = "https://www.dfam.org/releases/current/families/FamDB/"
# Place files where famdb.py expects them: <install_dir>/Libraries/famdb
# This script lives in <install_dir>/utils/, so go up one level.
_INSTALL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUT = os.path.join(_INSTALL_DIR, "Libraries", "famdb")

# Component display order and friendly names
COMPONENT_ORDER = ["root", "curated.consensus", "curated.hmm", "uncurated.consensus", "uncurated.hmm"]
COMPONENT_LABELS = {
    "root":                 "Root (taxonomy index)",
    "curated.consensus":    "Curated consensus sequences",
    "curated.hmm":          "Curated profile HMMs",
    "uncurated.consensus":  "Uncurated consensus sequences",
    "uncurated.hmm":        "Uncurated profile HMMs",
}


class _LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for attr, val in attrs:
                if attr == "href" and val and not val.startswith(("?", "/", "http")):
                    self.links.append(val.rstrip("/"))

    def error(self, message):
        pass


def _fetch_html(url):
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            return resp.read().decode("utf-8")
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach {url}: {e}") from e


def fetch_release_tag(release_root_url):
    """
    Look for a RELEASE_DFAM_X_Y tag file in the release root directory listing.
    Returns (tag, major, minor) or raises RuntimeError.
    """
    html = _fetch_html(release_root_url)

    parser = _LinkParser()
    parser.feed(html)

    release_pat = re.compile(r"^RELEASE_DFAM_(\d+)_(\d+)$")
    for link in parser.links:
        m = release_pat.match(link)
        if m:
            return link, int(m.group(1)), int(m.group(2))

    # Fallback: scan raw HTML text (handles non-linked entries)
    m = re.search(r"RELEASE_DFAM_(\d+)_(\d+)", html)
    if m:
        return m.group(0), int(m.group(1)), int(m.group(2))

    raise RuntimeError(f"No RELEASE_DFAM tag found at {release_root_url}")


def fetch_file_listing(base_url):
    """Return list of .h5.gz and .h5.gz.md5 filenames from the directory listing."""
    html = _fetch_html(base_url)
    parser = _LinkParser()
    parser.feed(html)
    return [f for f in parser.links if ".h5.gz" in f]


def parse_components(files):
    """
    Parse filenames into {component: [partition, ...]} and detect the release prefix.

    File patterns:
        {prefix}.{partition}.h5.gz              -- root
        {prefix}.{component}.{partition}.h5.gz  -- named component
    """
    gz_files = [f for f in files if f.endswith(".h5.gz")]
    root_pat = re.compile(r"^(dfam\w+?)\.(\d+)\.h5\.gz$")
    comp_pat = re.compile(r"^(dfam\w+?)\.((?:(?:un)?curated)\.(?:consensus|hmm))\.(\d+)\.h5\.gz$")

    prefix = None
    components = defaultdict(set)

    for f in gz_files:
        m = comp_pat.match(f)
        if m:
            prefix = m.group(1)
            components[m.group(2)].add(int(m.group(3)))
            continue
        m = root_pat.match(f)
        if m:
            prefix = m.group(1)
            components["root"].add(int(m.group(2)))

    return prefix, {k: sorted(v) for k, v in components.items()}


def _parse_range_selection(text, available_set):
    """Parse 'all', '0', '0-5', '0,2,5-10' against a set of available ints."""
    text = text.strip().lower()
    if text == "all":
        return sorted(available_set)

    result = set()
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            parts = token.split("-", 1)
            try:
                lo, hi = int(parts[0]), int(parts[1])
            except ValueError:
                continue
            result.update(p for p in available_set if lo <= p <= hi)
        else:
            try:
                p = int(token)
            except ValueError:
                continue
            if p in available_set:
                result.add(p)
            else:
                print(f"    Warning: partition {p} not available, skipping")
    return sorted(result)


def _prompt_components(components):
    """Interactive component selection; returns list of chosen component names."""
    ordered = [c for c in COMPONENT_ORDER if c in components]
    ordered += [c for c in sorted(components) if c not in ordered]

    print("\nAvailable components:")
    for i, comp in enumerate(ordered, 1):
        parts = components[comp]
        label = COMPONENT_LABELS.get(comp, comp)
        n = len(parts)
        part_range = f"partition {parts[0]}" if n == 1 else f"partitions 0-{parts[-1]} ({n} total)"
        print(f"  {i:2}. {label}")
        print(f"        [{comp}] -- {part_range}")

    print("\nEnter numbers to download (e.g. '1,3' or 'all'):")
    while True:
        raw = input("> ").strip().lower()
        if raw == "all":
            return ordered
        try:
            chosen = []
            for tok in raw.split(","):
                idx = int(tok.strip()) - 1
                if 0 <= idx < len(ordered):
                    chosen.append(ordered[idx])
            if chosen:
                return chosen
        except (ValueError, IndexError):
            pass
        print("Invalid selection, try again.")


def _prompt_partitions(comp, available):
    """Interactive partition selection for a single component."""
    label = COMPONENT_LABELS.get(comp, comp)

    if len(available) == 1:
        print(f"\n  {label}: only 1 partition ({available[0]}), selecting it.")
        return available

    print(f"\n  {label} has {len(available)} partitions (0-{available[-1]}).")
    print("  Enter partitions to download: 'all', '0', '0-5', '0,2,5-10', etc.")
    avail_set = set(available)
    while True:
        raw = input("  > ").strip()
        if not raw:
            continue
        result = _parse_range_selection(raw, avail_set)
        if result:
            return result
        print("  No valid partitions selected, try again.")


def _md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url, dest, label):
    """Download url -> dest with a simple progress line. Uses a .part temp file."""
    tmp = dest + ".part"
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            total = int(resp.headers.get("Content-Length", 0))
            done = 0
            with open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(65536)
                    if not chunk:
                        break
                    out.write(chunk)
                    done += len(chunk)
                    if total:
                        pct = done * 100 // total
                        mb = done / 1_048_576
                        print(f"\r    {label}: {pct:3d}%  {mb:.1f} MB", end="", flush=True)
                    else:
                        print(f"\r    {label}: {done / 1_048_576:.1f} MB", end="", flush=True)
        print()
        os.replace(tmp, dest)
    except Exception:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def _read_md5_file(path):
    """Read the first whitespace-delimited token from an md5 file."""
    with open(path) as f:
        return f.read().split()[0].lower()


def process_one(base_url, out_dir, prefix, component, partition, dry_run=False):
    """
    Download, validate, and decompress one partition file.
    Returns True on success (or already done), False on error.
    """
    stem = (
        f"{prefix}.{partition}.h5"
        if component == "root"
        else f"{prefix}.{component}.{partition}.h5"
    )
    gz_name  = stem + ".gz"
    md5_name = gz_name + ".md5"

    final = os.path.join(out_dir, stem)
    gz    = os.path.join(out_dir, gz_name)
    md5f  = os.path.join(out_dir, md5_name)

    tag = f"{component}/partition-{partition}"

    # Already fully decompressed?
    if os.path.exists(final):
        print(f"  {tag}: already present ({stem}), skipping.")
        return True

    if dry_run:
        print(f"  [dry-run] {gz_name}")
        return True

    # Fetch md5 sidecar
    if not os.path.exists(md5f):
        try:
            _download(base_url + md5_name, md5f, md5_name)
        except Exception as e:
            print(f"  ERROR fetching {md5_name}: {e}")
            return False

    try:
        expected = _read_md5_file(md5f)
    except Exception as e:
        print(f"  ERROR reading {md5f}: {e}")
        return False

    # Download gz if missing or corrupt
    need_dl = True
    if os.path.exists(gz):
        print(f"    Verifying existing {gz_name}...", end=" ", flush=True)
        if _md5sum(gz) == expected:
            print("OK")
            need_dl = False
        else:
            print("MISMATCH - re-downloading")
            os.remove(gz)

    if need_dl:
        try:
            _download(base_url + gz_name, gz, gz_name)
        except Exception as e:
            print(f"  ERROR downloading {gz_name}: {e}")
            return False

        actual = _md5sum(gz)
        if actual != expected:
            print(f"  ERROR: MD5 mismatch for {gz_name}")
            print(f"    expected: {expected}")
            print(f"    actual:   {actual}")
            os.remove(gz)
            return False
        print(f"    MD5 OK: {gz_name}")

    # Decompress
    print(f"    Decompressing {gz_name} ...", end=" ", flush=True)
    tmp_final = final + ".part"
    try:
        with gzip.open(gz, "rb") as f_in, open(tmp_final, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.replace(tmp_final, final)
        os.remove(gz)
        os.remove(md5f)
        print(f"done -> {stem}")
    except Exception as e:
        print(f"ERROR: {e}")
        for p in (tmp_final,):
            if os.path.exists(p):
                os.remove(p)
        return False

    return True


def main():
    ap = argparse.ArgumentParser(
        description="Download, verify, and decompress FamDB component files from Dfam."
    )
    ap.add_argument(
        "-o", "--output-dir",
        default=DEFAULT_OUT,
        metavar="DIR",
        help=f"Destination directory (default: {DEFAULT_OUT})",
    )
    ap.add_argument(
        "-u", "--url",
        default=DEFAULT_URL,
        metavar="URL",
        help=f"Base URL of FamDB release directory (default: {DEFAULT_URL})",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be downloaded without downloading anything",
    )
    args = ap.parse_args()

    base_url = args.url.rstrip("/") + "/"
    out_dir  = os.path.abspath(args.output_dir)

    # Navigate 2 directory levels up from the FamDB URL to find the release root
    release_root = base_url.rstrip("/").rsplit("/", 2)[0] + "/"
    try:
        tag, major, minor = fetch_release_tag(release_root)
    except RuntimeError as e:
        print(f"ERROR: Could not determine release version: {e}")
        sys.exit(1)

    print(f"Release: {tag}")

    MIN_MAJOR, MIN_MINOR = 4, 0
    if (major, minor) < (MIN_MAJOR, MIN_MINOR):
        print(f"ERROR: Release {tag} predates the minimum supported release "
              f"RELEASE_DFAM_{MIN_MAJOR}_{MIN_MINOR}.")
        print("Please report this problem to help@dfam.org")
        sys.exit(1)

    print(f"Fetching file listing from {base_url} ...")
    try:
        files = fetch_file_listing(base_url)
    except RuntimeError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    if not files:
        print("ERROR: No .h5.gz files found at the given URL.")
        sys.exit(1)

    prefix, components = parse_components(files)
    if not prefix:
        print("ERROR: Could not determine release prefix from filenames.")
        sys.exit(1)

    selected_components = _prompt_components(components)

    work = []
    for comp in selected_components:
        chosen = _prompt_partitions(comp, components[comp])
        for part in chosen:
            work.append((comp, part))

    if not work:
        print("Nothing selected.")
        sys.exit(0)

    print(f"\n{'[dry-run] ' if args.dry_run else ''}Downloading {len(work)} file(s) to: {out_dir}\n")
    if not args.dry_run:
        os.makedirs(out_dir, exist_ok=True)

    errors = []
    for i, (comp, part) in enumerate(work, 1):
        print(f"[{i}/{len(work)}] {comp} partition {part}")
        ok = process_one(base_url, out_dir, prefix, comp, part, dry_run=args.dry_run)
        if not ok:
            errors.append((comp, part))

    print("\n--- Summary ---")
    succeeded = len(work) - len(errors)
    print(f"  {succeeded}/{len(work)} completed successfully")
    if errors:
        print("  Failed:")
        for comp, part in errors:
            print(f"    {comp} partition {part}")
        sys.exit(1)


if __name__ == "__main__":
    main()
