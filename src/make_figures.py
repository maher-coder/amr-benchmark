"""
Generate 5 publication-quality figures for the CAMDA AMR paper.

Outputs to writeup/figures/:
  fig1_cross_source_drift.png — pct_R per source × combo heatmap
  fig2_cv_vs_test.png         — CV AUROC vs held-out AUROC with 95% CIs
  fig3_rule_vs_ml.png         — rule-only vs ML AUROC paired bars
  fig4_permutation_importance.png — top-5 features per combo
  fig5_scaling_curve.png      — AUROC vs N_train per combo (cleaner version)

Style: publication-friendly (clean, high-DPI, readable fonts).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import duckdb
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("Agg")
sys.stdout.reconfigure(line_buffering=True)

# Publication-friendly style
plt.rcParams.update({
    "font.size": 10,
    "font.family": "sans-serif",
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 12,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
})

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)

COMBOS = [
    ("Staphylococcus aureus", "methicillin", "S. aureus × methicillin"),
    ("Enterococcus faecium", "vancomycin", "E. faecium × vancomycin"),
    ("Neisseria gonorrhoeae", "ciprofloxacin", "N. gonorrhoeae × ciprofloxacin"),
    ("Escherichia coli", "ampicillin", "E. coli × ampicillin"),
    ("Salmonella enterica", "ampicillin", "Salmonella × ampicillin"),
    ("Salmonella enterica", "tetracycline", "Salmonella × tetracycline"),
]

COMBO_SHORT = {
    ("Staphylococcus aureus", "methicillin"): "Staph × meth",
    ("Enterococcus faecium", "vancomycin"): "E. faecium × vanco",
    ("Neisseria gonorrhoeae", "ciprofloxacin"): "N. gono × cipro",
    ("Escherichia coli", "ampicillin"): "E. coli × amp",
    ("Salmonella enterica", "ampicillin"): "Salmonella × amp",
    ("Salmonella enterica", "tetracycline"): "Salmonella × tet",
}


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — cross-source pct_R heatmap
# ──────────────────────────────────────────────────────────────────────────────
print("Fig 1 — cross-source drift heatmap", flush=True)
con = duckdb.connect(str(ROOT / "data" / "amr_portal" / "portal.duckdb"), read_only=True)

drift_rows = []
for sp_, drug, label in COMBOS:
    df = con.execute("""
        WITH expanded AS (
          SELECT assembly_ID, species, antibiotic_name, resistance_phenotype,
                 trim(s) AS src
          FROM phenotype, unnest(string_split(database, ';')) AS t(s)
          WHERE database IS NOT NULL
            AND species = ? AND antibiotic_name = ?
            AND resistance_phenotype IN ('resistant', 'susceptible')
        )
        SELECT src, COUNT(*) AS n,
               100.0 * SUM(CASE WHEN resistance_phenotype='resistant' THEN 1 ELSE 0 END) / COUNT(*) AS pct_r
        FROM expanded
        GROUP BY src
        HAVING COUNT(*) >= 50
        ORDER BY n DESC
    """, [sp_, drug]).fetchdf()
    for _, r in df.iterrows():
        drift_rows.append({
            "combo": COMBO_SHORT[(sp_, drug)],
            "source": r["src"],
            "pct_r": r["pct_r"],
            "n": r["n"],
        })
con.close()

drift_df = pd.DataFrame(drift_rows)
pivot = drift_df.pivot(index="source", columns="combo", values="pct_r")
# Keep same order as COMBOS
pivot = pivot[[COMBO_SHORT[(sp_, drug)] for sp_, drug, _ in COMBOS]]
# Sort sources by total sample count (rough)
source_order = ["PATRIC", "CABBAGE_PubMed_data", "NCBI_antibiogram", "NARMS",
                "NDARO", "COMPARE_ML", "pathogenwatch", "pubMLST", "microreact", "CDC"]
pivot = pivot.reindex([s for s in source_order if s in pivot.index])

fig, ax = plt.subplots(figsize=(10, 5.5))
mask = pivot.isna()
sns.heatmap(pivot, annot=True, fmt=".1f", cmap="RdYlBu_r", center=50,
            vmin=0, vmax=100, mask=mask, cbar_kws={"label": "% Resistant"},
            ax=ax, linewidths=0.4, linecolor="white",
            annot_kws={"fontsize": 8.5})
ax.set_xlabel("")
ax.set_ylabel("Source database")
ax.set_title("Resistance prevalence drift across cohort sources\n"
             "(N ≥ 50 per cell; grey = insufficient sample)")
ax.tick_params(axis="x", rotation=25)
plt.setp(ax.get_xticklabels(), ha="right")
ax.grid(False)
fig.tight_layout()
fig.savefig(OUT / "fig1_cross_source_drift.png")
plt.close(fig)
print(f"  → {OUT}/fig1_cross_source_drift.png", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — CV vs held-out test with 95% CIs
# ──────────────────────────────────────────────────────────────────────────────
print("Fig 2 — CV vs test AUROC", flush=True)
cv = pd.read_csv(ROOT / "results" / "cv_summary.csv")
test = pd.read_csv(ROOT / "results" / "held_out_test_metrics.csv")

# Parse 95% CI strings like "[0.9838, 0.9987]" -> (lo, hi)
def parse_ci(s):
    s = str(s).strip("[]")
    try:
        lo, hi = [float(x.strip()) for x in s.split(",")]
        return lo, hi
    except Exception:
        return np.nan, np.nan

test["auroc_lo"], test["auroc_hi"] = zip(*test["auroc_ci95"].map(parse_ci))
cv["auroc_lo"], cv["auroc_hi"] = zip(*cv["auroc_ci95"].map(parse_ci))

fig, ax = plt.subplots(figsize=(8.5, 6.5))

# Plot: EN model only for clarity
for i, (sp_, drug, label) in enumerate(COMBOS):
    cv_row = cv[(cv["species"] == sp_) & (cv["drug"] == drug) & (cv["model"] == "EN")]
    te_row = test[(test["species"] == sp_) & (test["drug"] == drug) & (test["model"] == "EN")]
    if len(cv_row) == 0 or len(te_row) == 0:
        continue
    cv_m = cv_row["auroc_mean"].iloc[0]
    cv_lo, cv_hi = cv_row["auroc_lo"].iloc[0], cv_row["auroc_hi"].iloc[0]
    te_m = te_row["auroc"].iloc[0]
    te_lo, te_hi = te_row["auroc_lo"].iloc[0], te_row["auroc_hi"].iloc[0]
    short = COMBO_SHORT[(sp_, drug)]
    ax.errorbar([cv_m], [te_m],
                xerr=[[cv_m - cv_lo], [cv_hi - cv_m]],
                yerr=[[te_m - te_lo], [te_hi - te_m]],
                fmt="o", markersize=9, capsize=4, linewidth=1.3,
                label=short, alpha=0.85)
    ax.annotate(short, (cv_m, te_m), xytext=(6, 6), textcoords="offset points", fontsize=8.5)

# Reference line y=x
xs = np.linspace(0.4, 1.0, 20)
ax.plot(xs, xs, "k--", linewidth=0.8, alpha=0.5, label="CV = Test")
# Shade ±0.05 band around y=x
ax.fill_between(xs, xs - 0.05, xs + 0.05, color="gray", alpha=0.1,
                label="±0.05 agreement band")

ax.set_xlabel("CV AUROC (train-pool GroupKFold, mean ± 95% CI)")
ax.set_ylabel("Held-out test AUROC (CABBAGE_PubMed, mean ± 95% CI)")
ax.set_title("Cross-validation vs independent held-out test generalization\n"
             "(ElasticNet, 6 pathogen-drug combinations)")
ax.set_xlim(0.55, 1.02)
ax.set_ylim(0.55, 1.02)
ax.legend(loc="lower right", fontsize=8, framealpha=0.85)
fig.tight_layout()
fig.savefig(OUT / "fig2_cv_vs_test.png")
plt.close(fig)
print(f"  → {OUT}/fig2_cv_vs_test.png", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — rule-only vs ML AUROC
# ──────────────────────────────────────────────────────────────────────────────
print("Fig 3 — rule-only vs ML AUROC", flush=True)
rule = pd.read_csv(ROOT / "results" / "rule_only_baseline.csv")

fig, ax = plt.subplots(figsize=(9, 5.2))
labels = []
rule_vals = []
ml_cv_vals = []
ml_te_vals = []
for sp_, drug, _ in COMBOS:
    rr = rule[(rule["species"] == sp_) & (rule["drug"] == drug)]
    cv_r = cv[(cv["species"] == sp_) & (cv["drug"] == drug) & (cv["model"] == "EN")]
    te_r = test[(test["species"] == sp_) & (test["drug"] == drug) & (test["model"] == "EN")]
    if len(rr) == 0:
        continue
    labels.append(COMBO_SHORT[(sp_, drug)])
    rule_vals.append(float(rr["train_auroc"].iloc[0]))
    ml_cv_vals.append(float(cv_r["auroc_mean"].iloc[0]))
    ml_te_vals.append(float(te_r["auroc"].iloc[0]))

x = np.arange(len(labels))
w = 0.26
ax.bar(x - w, rule_vals, w, label="Rule-only (canonical genes)", color="#888888")
ax.bar(x,       ml_cv_vals, w, label="ML ElasticNet — CV", color="#3b82f6")
ax.bar(x + w,   ml_te_vals, w, label="ML ElasticNet — held-out test", color="#1e40af")

for i, (r, c, t) in enumerate(zip(rule_vals, ml_cv_vals, ml_te_vals)):
    ax.text(i - w, r + 0.005, f"{r:.3f}", ha="center", fontsize=7.5)
    ax.text(i,     c + 0.005, f"{c:.3f}", ha="center", fontsize=7.5)
    ax.text(i + w, t + 0.005, f"{t:.3f}", ha="center", fontsize=7.5)

ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=20, ha="right")
ax.set_ylabel("AUROC")
ax.set_ylim(0.4, 1.05)
ax.set_title("ML value-add over canonical-gene rule baseline\n"
             "Grey bar = 'if gene present → resistant'; blue bars = ElasticNet")
ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
ax.axhline(y=0.5, color="red", linestyle=":", linewidth=0.6, alpha=0.5)
fig.tight_layout()
fig.savefig(OUT / "fig3_rule_vs_ml.png")
plt.close(fig)
print(f"  → {OUT}/fig3_rule_vs_ml.png", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — top-5 permutation importance per combo
# ──────────────────────────────────────────────────────────────────────────────
print("Fig 4 — permutation importance", flush=True)
perm = pd.read_csv(ROOT / "results" / "permutation_importance.csv")

fig, axes = plt.subplots(2, 3, figsize=(14, 7))
axes = axes.flatten()

for i, (sp_, drug, label) in enumerate(COMBOS):
    ax = axes[i]
    sub = perm[(perm["species"] == sp_) & (perm["drug"] == drug)]
    top = sub.sort_values("perm_drop_mean", ascending=False).head(5)
    top = top.sort_values("perm_drop_mean")  # ascending for horizontal bar display
    colors = plt.cm.viridis(np.linspace(0.2, 0.85, len(top)))
    ax.barh(top["feature"], top["perm_drop_mean"],
            xerr=top["perm_drop_std"], color=colors,
            edgecolor="#222222", linewidth=0.4, capsize=3)
    ax.set_xlabel("AUROC drop under permutation")
    ax.set_title(COMBO_SHORT[(sp_, drug)])
    ax.grid(axis="x", alpha=0.3)
    ax.tick_params(axis="y", labelsize=8.5)
    for j, v in enumerate(top["perm_drop_mean"].values):
        ax.text(v + 0.005 * max(0.01, top["perm_drop_mean"].max()),
                j, f"{v:.3f}", va="center", fontsize=7.5)

fig.suptitle("Top-5 features by permutation importance (honest metric vs LightGBM gain)",
             fontsize=12, y=1.00)
fig.tight_layout()
fig.savefig(OUT / "fig4_permutation_importance.png")
plt.close(fig)
print(f"  → {OUT}/fig4_permutation_importance.png", flush=True)


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 5 — scaling curve (cleaner, EN only)
# ──────────────────────────────────────────────────────────────────────────────
print("Fig 5 — scaling curve", flush=True)
scaling = pd.read_csv(ROOT / "results" / "scaling_curve.csv")

fig, ax = plt.subplots(figsize=(9, 5.5))
colors = plt.cm.tab10(np.linspace(0, 1, 10))
for i, (sp_, drug, label) in enumerate(COMBOS):
    sub = scaling[(scaling["species"] == sp_) & (scaling["drug"] == drug) & (scaling["model"] == "EN")]
    sub = sub.sort_values("n_actual")
    ax.plot(sub["n_actual"], sub["auroc"], marker="o", linewidth=1.8,
            markersize=6, label=COMBO_SHORT[(sp_, drug)], color=colors[i])

ax.set_xscale("log")
ax.set_xlabel("Training set size (N)")
ax.set_ylabel("AUROC (GroupKFold-by-source, mean)")
ax.set_title("Scaling curves — AUROC vs training data volume\n"
             "(ElasticNet, train-pool only; flat plateau = feature-limited, not data-limited)")
ax.set_ylim(0.4, 1.02)
ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
fig.tight_layout()
fig.savefig(OUT / "fig5_scaling_curve.png")
plt.close(fig)
print(f"  → {OUT}/fig5_scaling_curve.png", flush=True)

print("\nAll 5 figures saved to writeup/figures/", flush=True)
