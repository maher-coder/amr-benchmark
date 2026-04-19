"""
Fetch the CABBAGE December 2025 release from the EBI AMR portal.

Public FTP, no DUA or credentialing required. Downloads the phenotype /
genotype / merged parquet files plus the full DuckDB portal (167 MB).

After fetch:
    data/amr_portal/portal.duckdb         (167 MB — fully queryable)
    data/amr_portal/phenotype.parquet      (8.3 MB)
    data/amr_portal/genotype.parquet       (22 MB)
    data/amr_portal/phenotype_genotype_merged.parquet  (4.3 MB)
    data/amr_portal/manifest.yaml          (tiny)
    data/amr_portal/md5sums                (tiny)

Usage:
    python src/fetch_data.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/amr_portal/releases/2025-12"
FILES = [
    "portal.duckdb",
    "phenotype.parquet",
    "genotype.parquet",
    "phenotype_genotype_merged.parquet",
    "manifest.yaml",
    "md5sums",
]


def main():
    data_dir = Path(__file__).resolve().parent.parent / "data" / "amr_portal"
    data_dir.mkdir(parents=True, exist_ok=True)
    for fname in FILES:
        target = data_dir / fname
        if target.exists():
            print(f"✓ {fname} exists ({target.stat().st_size / 1e6:.1f} MB)", flush=True)
            continue
        url = f"{BASE_URL}/{fname}"
        print(f"→ fetching {url}", flush=True)
        subprocess.run(["wget", "-q", "-O", str(target), url], check=True)
        print(f"✓ {fname} downloaded ({target.stat().st_size / 1e6:.1f} MB)", flush=True)

    # Verify md5
    print("\n→ verifying md5 checksums...", flush=True)
    result = subprocess.run(
        ["md5sum", "-c", "md5sums"], cwd=data_dir, capture_output=True, text=True
    )
    for line in result.stdout.splitlines() + result.stderr.splitlines():
        if "OK" in line or "FAIL" in line:
            print(f"   {line}", flush=True)
    print("\nDownload complete. Next step: python src/prepare.py", flush=True)


if __name__ == "__main__":
    main()
