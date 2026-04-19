"""
Train ElasticNet and LightGBM per combo with GroupKFold-by-source CV,
reserve CABBAGE_PubMed_data as held-out test cohort, evaluate once at the end.

Outputs:
    results/cv_summary.csv            — per-combo CV metrics with bootstrap CI
    results/held_out_test_metrics.csv — final test metrics (single-shot)
    results/rule_only_baseline.csv    — rule-only comparison

Usage:
    python src/train.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import GroupKFold, StratifiedKFold

sys.stdout.reconfigure(line_buffering=True)
warnings.filterwarnings("ignore")

import lightgbm as lgb

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

SEED = 42
N_FOLDS = 5
HELD_OUT_SOURCE = "CABBAGE_PubMed_data"

COMBOS = [
    ("Staphylococcus aureus", "methicillin"),
    ("Enterococcus faecium", "vancomycin"),
    ("Neisseria gonorrhoeae", "ciprofloxacin"),
    ("Escherichia coli", "ampicillin"),
    ("Salmonella enterica", "ampicillin"),
    ("Salmonella enterica", "tetracycline"),
]

CANONICAL_RULES = {
    ("Staphylococcus aureus", "methicillin"): ["mecA"],
    ("Enterococcus faecium", "vancomycin"): ["vanA", "vanB", "vanX-A", "vanR-A", "vanH-A", "vanY-A"],
    ("Neisseria gonorrhoeae", "ciprofloxacin"): ["gyrA_S91F", "gyrA_S83L", "parC_S87R", "parC_D86N"],
    ("Escherichia coli", "ampicillin"): ["blaTEM-1", "blaCMY-2", "blaCTX-M-15", "blaCTX-M-14", "blaCTX-M-1", "blaCTX-M-27", "blaCTX-M-55", "blaSHV-1", "blaOXA-48", "blaHER-3"],
    ("Salmonella enterica", "ampicillin"): ["blaTEM-1", "blaCMY-2", "blaCTX-M-65", "blaCTX-M-1", "blaCTX-M-14", "blaTEM-135", "blaHER-3"],
    ("Salmonella enterica", "tetracycline"): ["tet(A)", "tet(B)", "tet(C)", "tet(G)", "tet(D)"],
}


def load():
    X = sp.load_npz(DATA / "feature_matrix.npz").astype(np.float32)
    rows = pd.read_csv(DATA / "feature_matrix_rows.csv")["assembly_ID"].tolist()
    cols = pd.read_csv(DATA / "feature_matrix_cols.csv")["amr_element_symbol"].tolist()
    pheno = pd.read_parquet(DATA / "phenotype_clean.parquet")
    pheno["primary_source"] = (
        pheno["database"].fillna("").str.split(";").str[0].str.strip().replace("", "unknown")
    )
    return X, rows, cols, pheno


def make_models():
    en = LogisticRegression(penalty="l2", C=1.0, class_weight="balanced",
                            max_iter=500, solver="liblinear", random_state=SEED)
    lgbm = lgb.LGBMClassifier(n_estimators=300, learning_rate=0.05, num_leaves=31,
                              min_child_samples=10, reg_alpha=0.1, reg_lambda=0.1,
                              subsample=0.8, colsample_bytree=0.8, objective="binary",
                              is_unbalance=True, random_state=SEED, n_jobs=-1, verbose=-1)
    return {"EN": en, "LGB": lgbm}


def bootstrap_ci(values, n_boot=1000, seed=SEED):
    rng = np.random.default_rng(seed)
    boots = [rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_boot)]
    return np.mean(values), *np.quantile(boots, [0.025, 0.975])


def bootstrap_metric(y_true, y_score, metric_fn, n_boot=1000, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(y_true), size=len(y_true))
        try:
            vals.append(metric_fn(y_true[idx], y_score[idx]))
        except ValueError:
            continue
    return np.quantile(vals, [0.025, 0.975])


def cv_combo(X, rows, pheno_train, sp_, drug):
    asm_to_idx = {a: i for i, a in enumerate(rows)}
    sub = pheno_train[(pheno_train["species"] == sp_) & (pheno_train["antibiotic_name"] == drug)].copy()
    sub["y_bin"] = (sub["y"] == "R").astype(int)
    sub = sub[sub["assembly_ID"].isin(asm_to_idx)]
    if len(sub) < 100: return None
    idxs = sub["assembly_ID"].map(asm_to_idx).values
    Xc = X[idxs]
    y = sub["y_bin"].values
    groups = sub["primary_source"].values
    n_groups = len(np.unique(groups))
    if n_groups >= 2:
        splitter = GroupKFold(n_splits=min(N_FOLDS, n_groups))
        splits = list(splitter.split(Xc, y, groups=groups))
    else:
        splitter = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
        splits = list(splitter.split(Xc, y))

    per_fold = []
    for fold, (tr, va) in enumerate(splits):
        if len(np.unique(y[tr])) < 2: continue
        models = make_models()
        fold_row = {"fold": fold, "n_va": len(va)}
        for name, m in models.items():
            m.fit(Xc[tr], y[tr])
            pp = m.predict_proba(Xc[va])[:, 1]
            p = m.predict(Xc[va])
            try:
                fold_row[f"{name}_auroc"] = roc_auc_score(y[va], pp)
            except ValueError:
                fold_row[f"{name}_auroc"] = np.nan
            fold_row[f"{name}_mcc"] = matthews_corrcoef(y[va], p)
            fold_row[f"{name}_bacc"] = balanced_accuracy_score(y[va], p)
        per_fold.append(fold_row)
    return per_fold


def evaluate_test(X, rows, pheno_train, pheno_test, sp_, drug, model_name):
    asm_to_idx = {a: i for i, a in enumerate(rows)}
    tr = pheno_train[(pheno_train["species"] == sp_) & (pheno_train["antibiotic_name"] == drug)].copy()
    te = pheno_test[(pheno_test["species"] == sp_) & (pheno_test["antibiotic_name"] == drug)].copy()
    tr["y_bin"] = (tr["y"] == "R").astype(int)
    te["y_bin"] = (te["y"] == "R").astype(int)
    tr = tr[tr["assembly_ID"].isin(asm_to_idx)]
    te = te[te["assembly_ID"].isin(asm_to_idx)]
    if len(tr) < 50 or len(te) < 20: return None

    Xtr, Xte = X[tr["assembly_ID"].map(asm_to_idx).values], X[te["assembly_ID"].map(asm_to_idx).values]
    ytr, yte = tr["y_bin"].values, te["y_bin"].values
    m = make_models()[model_name]
    m.fit(Xtr, ytr)
    pp, p = m.predict_proba(Xte)[:, 1], m.predict(Xte)

    try:
        auroc = roc_auc_score(yte, pp)
        auroc_ci = bootstrap_metric(yte, pp, roc_auc_score)
    except ValueError:
        auroc, auroc_ci = np.nan, (np.nan, np.nan)

    return {
        "species": sp_, "drug": drug, "model": model_name,
        "n_train": len(tr), "n_test": len(te), "test_r_pct": round(100 * yte.mean(), 1),
        "auroc": round(auroc, 4), "auroc_ci_lo": round(auroc_ci[0], 4), "auroc_ci_hi": round(auroc_ci[1], 4),
        "mcc": round(matthews_corrcoef(yte, p), 4),
        "bacc": round(balanced_accuracy_score(yte, p), 4),
        "acc": round(accuracy_score(yte, p), 4),
    }


def rule_baseline(X, rows, cols, pheno_train, sp_, drug):
    genes = CANONICAL_RULES[(sp_, drug)]
    gene_idx = [cols.index(g) for g in genes if g in cols]
    asm_to_idx = {a: i for i, a in enumerate(rows)}
    sub = pheno_train[(pheno_train["species"] == sp_) & (pheno_train["antibiotic_name"] == drug)].copy()
    sub["y_bin"] = (sub["y"] == "R").astype(int)
    sub = sub[sub["assembly_ID"].isin(asm_to_idx)]
    if len(sub) < 50: return None
    idxs = sub["assembly_ID"].map(asm_to_idx).values
    X_sub = X[idxs]
    y = sub["y_bin"].values
    score = X_sub[:, gene_idx].sum(axis=1).A1.astype(float)
    pred = (score > 0).astype(int)
    try:
        auroc = roc_auc_score(y, score)
    except ValueError:
        auroc = np.nan
    return {
        "species": sp_, "drug": drug, "n": len(sub),
        "auroc": round(auroc, 4),
        "bacc": round(balanced_accuracy_score(y, pred), 4),
        "mcc": round(matthews_corrcoef(y, pred), 4),
        "accuracy": round(accuracy_score(y, pred), 4),
    }


def main():
    print("=" * 70 + "\nTrain + evaluate\n" + "=" * 70, flush=True)
    X, rows, cols, pheno = load()
    pheno_train = pheno[pheno["primary_source"] != HELD_OUT_SOURCE].copy()
    pheno_test = pheno[pheno["primary_source"] == HELD_OUT_SOURCE].copy()
    print(f"train pool: {len(pheno_train):,}   held-out: {len(pheno_test):,}", flush=True)

    # CV
    cv_rows = []
    for sp_, drug in COMBOS:
        print(f"\n>> CV: {sp_} × {drug}", flush=True)
        folds = cv_combo(X, rows, pheno_train, sp_, drug)
        if folds is None: continue
        for f in folds:
            f["species"], f["drug"] = sp_, drug
            cv_rows.append(f)
    cv_df = pd.DataFrame(cv_rows)
    cv_summary = []
    for (sp_, drug), g in cv_df.groupby(["species", "drug"]):
        for model in ["EN", "LGB"]:
            aur = g[f"{model}_auroc"].dropna().values
            mean, lo, hi = bootstrap_ci(aur) if len(aur) else (np.nan, np.nan, np.nan)
            cv_summary.append({
                "species": sp_, "drug": drug, "model": model,
                "n_folds": len(g),
                "auroc_mean": round(mean, 4),
                "auroc_std": round(float(np.std(aur)), 4) if len(aur) else np.nan,
                "auroc_ci_lo": round(lo, 4), "auroc_ci_hi": round(hi, 4),
                "mcc_mean": round(float(g[f"{model}_mcc"].mean()), 4),
                "bacc_mean": round(float(g[f"{model}_bacc"].mean()), 4),
            })
    pd.DataFrame(cv_summary).to_csv(RESULTS / "cv_summary.csv", index=False)
    print(f"\n→ results/cv_summary.csv saved", flush=True)

    # Held-out single-shot
    print("\n" + "=" * 70 + "\nHELD-OUT SINGLE-SHOT EVALUATION\n" + "=" * 70, flush=True)
    test_rows = []
    for sp_, drug in COMBOS:
        for model in ["EN", "LGB"]:
            res = evaluate_test(X, rows, pheno_train, pheno_test, sp_, drug, model)
            if res:
                test_rows.append(res)
                print(f"{sp_} × {drug} [{model}]: AUROC {res['auroc']} [{res['auroc_ci_lo']}, {res['auroc_ci_hi']}] "
                      f"MCC {res['mcc']} Bacc {res['bacc']}", flush=True)
    pd.DataFrame(test_rows).to_csv(RESULTS / "held_out_test_metrics.csv", index=False)
    print(f"→ results/held_out_test_metrics.csv saved", flush=True)

    # Rule-only baseline
    print("\n" + "=" * 70 + "\nRULE-ONLY BASELINE\n" + "=" * 70, flush=True)
    rule_rows = []
    for sp_, drug in COMBOS:
        res = rule_baseline(X, rows, cols, pheno_train, sp_, drug)
        if res:
            rule_rows.append(res)
            print(f"{sp_} × {drug}: AUROC {res['auroc']}  MCC {res['mcc']}", flush=True)
    pd.DataFrame(rule_rows).to_csv(RESULTS / "rule_only_baseline.csv", index=False)
    print(f"→ results/rule_only_baseline.csv saved", flush=True)
    print("\nAll results in results/.", flush=True)


if __name__ == "__main__":
    main()
