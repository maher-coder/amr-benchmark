"""
Prepare the CABBAGE data for training.

Applies six filters documented in the paper (Methods §2.2):
  1. assembly_ID NOT NULL  (removes 156k phenotype-only rows)
  2. phenotype = COALESCE(Updated_EUCAST, Updated_CLSI, legacy) ∈ {R, S}
  3. collection_year ≥ 1945 OR NULL  (removes 2 200 rows dated 1922)
  4. filter biologically impossible combinations
     (vanco × Gram-neg; colistin × Gram+; metronidazole × aerobes)
  5. dedupe (assembly_ID, antibiotic_name) via majority vote (drop ties)
  6. build binary presence/absence matrix over amr_element_symbol, min_n=50

Outputs:
    data/phenotype_clean.parquet        — filtered labels
    data/feature_matrix.npz             — sparse CSR binary
    data/feature_matrix_rows.csv        — assembly_ID order
    data/feature_matrix_cols.csv        — amr_element_symbol order

Usage:
    python src/prepare.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.stdout.reconfigure(line_buffering=True)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DUCKDB_PATH = DATA / "amr_portal" / "portal.duckdb"
MIN_N_FEATURE = 50


def main():
    print("=" * 70, flush=True)
    print("CABBAGE data preparation", flush=True)
    print("=" * 70, flush=True)

    con = duckdb.connect(str(DUCKDB_PATH), read_only=True)

    # 1-5: Apply filters + dedup via SQL + Python
    print("\n[1/2] Applying 6 filters + majority-vote dedup...", flush=True)
    raw = con.execute("""
        WITH coalesced AS (
          SELECT
            assembly_ID, species, antibiotic_name, database, collection_year,
            COALESCE(
              CASE Updated_phenotype_EUCAST WHEN 'resistant' THEN 'R'
                                            WHEN 'susceptible' THEN 'S' END,
              CASE Updated_phenotype_CLSI   WHEN 'resistant' THEN 'R'
                                            WHEN 'susceptible' THEN 'S' END,
              CASE resistance_phenotype     WHEN 'resistant' THEN 'R'
                                            WHEN 'susceptible' THEN 'S' END
            ) AS y
          FROM phenotype
          WHERE assembly_ID IS NOT NULL
            AND (collection_year >= 1945 OR collection_year IS NULL)
            AND NOT (antibiotic_name = 'vancomycin'
                     AND species IN ('Escherichia coli','Klebsiella pneumoniae',
                                     'Pseudomonas aeruginosa'))
            AND NOT (antibiotic_name = 'colistin' AND species = 'Staphylococcus aureus')
            AND NOT (antibiotic_name = 'metronidazole'
                     AND species IN ('Escherichia coli','Klebsiella pneumoniae',
                                     'Pseudomonas aeruginosa','Staphylococcus aureus'))
        )
        SELECT assembly_ID, species, antibiotic_name, database, y
        FROM coalesced
        WHERE y IN ('R', 'S')
    """).fetchdf()
    print(f"  rows before dedup: {len(raw):,}", flush=True)

    # Majority-vote dedup of (assembly_ID, antibiotic_name)
    pheno = (
        raw.groupby(["assembly_ID", "antibiotic_name"])
        .agg(y=("y", lambda s: s.mode().iloc[0] if len(s.mode()) == 1 else np.nan),
             species=("species", "first"),
             database=("database", "first"),
             n_rows=("y", "size"))
        .reset_index()
        .dropna(subset=["y"])
    )
    print(f"  rows after dedup: {len(pheno):,}", flush=True)
    print(f"  unique assemblies: {pheno['assembly_ID'].nunique():,}", flush=True)

    pheno.to_parquet(DATA / "phenotype_clean.parquet", index=False)
    print(f"  → saved data/phenotype_clean.parquet", flush=True)

    # 6: Build feature matrix
    print(f"\n[2/2] Building binary feature matrix (min_n={MIN_N_FEATURE})...", flush=True)
    labeled_asm = pheno["assembly_ID"].unique().tolist()
    con.execute("""
        CREATE OR REPLACE TEMP TABLE labeled_asm AS
        SELECT DISTINCT unnest(?::VARCHAR[]) AS assembly_ID
    """, [labeled_asm])

    top_elements = con.execute(f"""
        SELECT amr_element_symbol, COUNT(DISTINCT assembly_ID) AS n_asm
        FROM genotype
        WHERE amr_element_symbol IS NOT NULL AND amr_element_symbol != ''
          AND assembly_ID IN (SELECT assembly_ID FROM labeled_asm)
        GROUP BY amr_element_symbol
        HAVING n_asm >= {MIN_N_FEATURE}
        ORDER BY n_asm DESC
    """).fetchdf()
    print(f"  features with N≥{MIN_N_FEATURE}: {len(top_elements):,}", flush=True)

    gene_pres = con.execute(f"""
        SELECT DISTINCT assembly_ID, amr_element_symbol
        FROM genotype
        WHERE assembly_ID IN (SELECT assembly_ID FROM labeled_asm)
          AND amr_element_symbol IN (
            SELECT amr_element_symbol FROM (
              SELECT amr_element_symbol, COUNT(DISTINCT assembly_ID) AS n
              FROM genotype
              WHERE amr_element_symbol IS NOT NULL AND amr_element_symbol != ''
                AND assembly_ID IN (SELECT assembly_ID FROM labeled_asm)
              GROUP BY amr_element_symbol
              HAVING n >= {MIN_N_FEATURE}
            )
          )
    """).fetchdf()

    asm_ids = sorted(labeled_asm)
    elem_ids = sorted(top_elements["amr_element_symbol"].unique())
    asm_to_idx = {a: i for i, a in enumerate(asm_ids)}
    elem_to_idx = {e: i for i, e in enumerate(elem_ids)}

    rows = gene_pres["assembly_ID"].map(asm_to_idx).values
    cols = gene_pres["amr_element_symbol"].map(elem_to_idx).values
    data_v = np.ones(len(gene_pres), dtype=np.uint8)
    X = sp.csr_matrix((data_v, (rows, cols)), shape=(len(asm_ids), len(elem_ids)))
    print(f"  matrix shape: {X.shape}  density: {X.nnz / (X.shape[0] * X.shape[1]):.4f}", flush=True)

    sp.save_npz(DATA / "feature_matrix.npz", X)
    pd.Series(asm_ids, name="assembly_ID").to_csv(DATA / "feature_matrix_rows.csv", index=False)
    pd.Series(elem_ids, name="amr_element_symbol").to_csv(DATA / "feature_matrix_cols.csv", index=False)
    print(f"  → saved data/feature_matrix.npz + row/col CSVs", flush=True)

    con.close()
    print("\nPrepare complete. Next step: python src/train.py", flush=True)


if __name__ == "__main__":
    main()
