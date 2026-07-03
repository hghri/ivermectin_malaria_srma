#!/usr/bin/env python3
"""
Python reproducibility pipeline for the ivermectin-malaria systematic review.

Purpose
-------
1. Validate the final extraction workbook.
2. Reproduce outcome-specific REML meta-analysis summaries.
3. Produce publication-ready summary and individual forest plots.
4. Export clean CSV tables for manuscript and appendix.

This Python implementation is intended as an independent check and plotting
utility. The primary statistical analysis for publication should be run using
the R scripts, which use metafor and clubSandwich.

Expected input
--------------
data/ivermectin_FINAL.xlsx with a sheet named Poolable_post_MDA.

Key required columns
--------------------
Author_Year, Country, Domain, Outcome_measure, Timepoint_clean, Comparison,
yi, sei, vi, poolable_primary.

Outputs
-------
outputs/tables/python_reml_outcome_summary.csv
outputs/figures/Figure3_primary_outcome_REML_forest.png/pdf
outputs/figures/forest_individual/*.png/pdf
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "ivermectin_FINAL.xlsx"
OUT_TABLES = ROOT / "outputs" / "tables"
OUT_FIGS = ROOT / "outputs" / "figures"
OUT_INDIVIDUAL = OUT_FIGS / "forest_individual"

OUT_TABLES.mkdir(parents=True, exist_ok=True)
OUT_INDIVIDUAL.mkdir(parents=True, exist_ok=True)

REQUIRED_COLUMNS = [
    "Author_Year", "Country", "Domain", "Outcome_measure", "Timepoint_clean",
    "Comparison", "yi", "sei", "vi", "poolable_primary"
]

OUTCOME_ORDER = [
    "Malaria incidence/infection",
    "Entomological inoculation rate",
    "Mosquito mortality",
    "Mosquito survival",
    "Vector density",
    "Sporozoite rate",
    "Parity rate",
    "Human biting rate",
    "Blood-feeding proportion",
    "Mosquito fecundity",
    "Any adverse events",
    "Serious adverse events/deaths",
    "Other outcome",
]

DOMAIN_ORDER = {"Malaria epidemiology": 1, "Entomological": 2, "Safety": 3}


def clean_filename(text: str) -> str:
    text = text.lower().strip()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def load_poolable_data(path: Path = DATA) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input workbook not found: {path}")
    df = pd.read_excel(path, sheet_name="Poolable_post_MDA")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in Poolable_post_MDA: {missing}")

    # Keep primary post-MDA poolable rows.
    # The workbook stores booleans as bool or strings depending on Excel export.
    pool = df["poolable_primary"].astype(str).str.lower().isin(["true", "1", "yes"])
    df = df.loc[pool].copy()

    # Ensure numeric effect and variance columns.
    for c in ["yi", "sei", "vi"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["yi", "vi"])
    df = df[df["vi"] >= 0].copy()
    df["sei"] = np.sqrt(df["vi"])

    # Clean display labels.
    for c in ["Author_Year", "Country", "Domain", "Outcome_measure", "Timepoint_clean", "Comparison"]:
        df[c] = df[c].fillna("").astype(str)
    df["study_label"] = (
        df["Author_Year"] + " | " + df["Country"] + " | " + df["Timepoint_clean"]
    ).str.replace(r"\s+", " ", regex=True).str.strip(" |")
    return df


def reml_intercept(y: np.ndarray, vi: np.ndarray) -> Dict[str, float]:
    """Fit random-effects intercept-only model using REML.

    This is a compact independent implementation used for validation and plotting.
    For final manuscript analysis, use R/metafor in R/02_reml_meta_analysis.R.
    """
    y = np.asarray(y, dtype=float)
    vi = np.asarray(vi, dtype=float)
    ok = np.isfinite(y) & np.isfinite(vi) & (vi >= 0)
    y, vi = y[ok], vi[ok]
    k = len(y)
    if k == 0:
        raise ValueError("No valid effects")
    if k == 1:
        mu = float(y[0])
        se = float(math.sqrt(max(vi[0], 0)))
        return {
            "k_effects": 1, "mu": mu, "se": se, "ci_l": mu - 1.96 * se,
            "ci_u": mu + 1.96 * se, "tau2": 0.0, "I2": 0.0,
            "pred_l": np.nan, "pred_u": np.nan, "Q": np.nan, "p_Q": np.nan
        }

    def nll(tau2: float) -> float:
        V = vi + tau2
        if np.any(V <= 0):
            return np.inf
        w = 1.0 / V
        sw = np.sum(w)
        mu = np.sum(w * y) / sw
        resid = y - mu
        # REML negative log-likelihood for intercept-only model.
        return 0.5 * (np.sum(np.log(V)) + np.log(sw) + np.sum(w * resid ** 2))

    upper = max(1.0, float(np.var(y) * 10 + np.max(vi) * 10))
    opt = minimize_scalar(nll, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-10})
    tau2 = max(0.0, float(opt.x))
    V = vi + tau2
    w = 1.0 / V
    mu = float(np.sum(w * y) / np.sum(w))
    se = float(math.sqrt(1.0 / np.sum(w)))

    # Cochran's Q using fixed-effect weights for reporting.
    wf = 1.0 / vi
    muf = float(np.sum(wf * y) / np.sum(wf))
    Q = float(np.sum(wf * (y - muf) ** 2))
    # I2 using metafor-style typical within-study variance approximation.
    sw = np.sum(wf)
    sw2 = np.sum(wf ** 2)
    typical_v = float((k - 1) * sw / (sw ** 2 - sw2)) if sw ** 2 > sw2 else float(np.mean(vi))
    I2 = max(0.0, 100.0 * tau2 / (tau2 + typical_v)) if tau2 + typical_v > 0 else 0.0

    return {
        "k_effects": k,
        "mu": mu,
        "se": se,
        "ci_l": mu - 1.96 * se,
        "ci_u": mu + 1.96 * se,
        "tau2": tau2,
        "I2": I2,
        "pred_l": mu - 1.96 * math.sqrt(tau2 + se ** 2),
        "pred_u": mu + 1.96 * math.sqrt(tau2 + se ** 2),
        "Q": Q,
        "p_Q": np.nan,
    }


def summarise_by_outcome(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for outcome, g in df.groupby("Outcome_measure", sort=False):
        fit = reml_intercept(g["yi"].values, g["vi"].values)
        domain = g["Domain"].mode().iloc[0] if not g["Domain"].mode().empty else ""
        rows.append({
            "Domain": domain,
            "Outcome": outcome,
            "Studies": g["Study_ID"].nunique() if "Study_ID" in g else g["Author_Year"].nunique(),
            "Effects": len(g),
            "Ratio": math.exp(fit["mu"]),
            "CI_lower": math.exp(fit["ci_l"]),
            "CI_upper": math.exp(fit["ci_u"]),
            "PI_lower": math.exp(fit["pred_l"]) if np.isfinite(fit["pred_l"]) else np.nan,
            "PI_upper": math.exp(fit["pred_u"]) if np.isfinite(fit["pred_u"]) else np.nan,
            "tau2": fit["tau2"],
            "I2": fit["I2"],
            "mu_log": fit["mu"],
            "se_log": fit["se"],
        })
    out = pd.DataFrame(rows)
    out["_domain_order"] = out["Domain"].map(DOMAIN_ORDER).fillna(99)
    out["_outcome_order"] = out["Outcome"].apply(lambda x: OUTCOME_ORDER.index(x) if x in OUTCOME_ORDER else 99)
    out = out.sort_values(["_domain_order", "_outcome_order", "Outcome"]).drop(columns=["_domain_order", "_outcome_order"])
    return out


def ratio_text(row: pd.Series) -> str:
    return f'{row["Ratio"]:.2f} ({row["CI_lower"]:.2f}–{row["CI_upper"]:.2f})'


def pi_text(row: pd.Series) -> str:
    if not np.isfinite(row["PI_lower"]):
        return "Not estimable"
    return f'{row["PI_lower"]:.2f}–{row["PI_upper"]:.2f}'


def plot_summary_forest(summary: pd.DataFrame) -> None:
    plot_df = summary.copy()
    plot_df["label"] = plot_df["Outcome"]
    plot_df = plot_df.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_df))

    fig_h = max(7, 0.55 * len(plot_df) + 2)
    fig, ax = plt.subplots(figsize=(9.5, fig_h))

    x = plot_df["Ratio"].values
    xl = plot_df["CI_lower"].values
    xu = plot_df["CI_upper"].values
    ax.errorbar(x, y, xerr=[x - xl, xu - x], fmt="o", capsize=3, markersize=5)

    # Add prediction interval as a thicker, lighter horizontal segment where available.
    for i, row in plot_df.iterrows():
        if np.isfinite(row["PI_lower"]):
            ax.hlines(y=i - 0.18, xmin=row["PI_lower"], xmax=row["PI_upper"], linewidth=4, alpha=0.25)

    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["label"], fontsize=9)
    ax.set_xlabel("Ratio estimate (log scale); values <1 favour ivermectin for malaria/entomological outcomes")
    ax.set_title("Figure 3. Primary outcome-specific REML random-effects meta-analysis")
    ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.7)

    # Right-side text columns.
    xlim = ax.get_xlim()
    right_x = xlim[1] * 1.15
    ax.set_xlim(xlim[0], xlim[1] * 2.2)
    for i, row in plot_df.iterrows():
        ax.text(right_x, i, f'{ratio_text(row)}; PI {pi_text(row)}', va="center", fontsize=8)
    ax.text(right_x, len(plot_df) + 0.3, "REML ratio (95% CI); prediction interval", fontsize=8, fontweight="bold")

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(OUT_FIGS / f"Figure3_primary_outcome_REML_forest.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_individual_forest(df: pd.DataFrame, outcome: str) -> Optional[Tuple[Path, Path]]:
    g = df[df["Outcome_measure"] == outcome].copy()
    if g.empty:
        return None
    fit = reml_intercept(g["yi"].values, g["vi"].values)
    g["ratio"] = np.exp(g["yi"])
    g["ci_l"] = np.exp(g["yi"] - 1.96 * np.sqrt(g["vi"]))
    g["ci_u"] = np.exp(g["yi"] + 1.96 * np.sqrt(g["vi"]))
    g = g.sort_values(["Author_Year", "Country", "Timepoint_clean", "Comparison"]).reset_index(drop=True)

    labels = g["study_label"].str.slice(0, 80).tolist()
    n = len(g)
    fig_h = max(5, min(18, 0.32 * n + 3.5))
    fig, ax = plt.subplots(figsize=(10, fig_h))
    y = np.arange(n, 0, -1)

    x = g["ratio"].values
    xl = g["ci_l"].values
    xu = g["ci_u"].values
    ax.errorbar(x, y, xerr=[x - xl, xu - x], fmt="o", capsize=2, markersize=3.5, linewidth=0.8)

    # Pooled estimate as diamond-like errorbar at y=0.
    pooled = math.exp(fit["mu"])
    ci_l = math.exp(fit["ci_l"])
    ci_u = math.exp(fit["ci_u"])
    ax.errorbar([pooled], [0], xerr=[[pooled - ci_l], [ci_u - pooled]], fmt="D", markersize=7, capsize=4, linewidth=1.5)

    # Prediction interval.
    if np.isfinite(fit["pred_l"]):
        ax.hlines(-0.35, xmin=math.exp(fit["pred_l"]), xmax=math.exp(fit["pred_u"]), linewidth=4, alpha=0.25)

    ax.axvline(1.0, linestyle="--", linewidth=1)
    ax.set_xscale("log")
    ax.set_yticks(list(y) + [0])
    ax.set_yticklabels(labels + ["Pooled REML"], fontsize=7)
    ax.set_xlabel("Ratio estimate (log scale)")
    title = f"{outcome}: individual study estimates and REML pooled effect"
    ax.set_title(title, fontsize=11)
    ax.grid(axis="x", linestyle=":", linewidth=0.5, alpha=0.7)

    note = f"Pooled {pooled:.2f} ({ci_l:.2f}–{ci_u:.2f}); tau²={fit['tau2']:.3f}; I²={fit['I2']:.1f}%"
    if np.isfinite(fit["pred_l"]):
        note += f"; PI {math.exp(fit['pred_l']):.2f}–{math.exp(fit['pred_u']):.2f}"
    ax.text(0.01, -0.12, note, transform=ax.transAxes, fontsize=8, va="top")

    fig.tight_layout()
    base = clean_filename(outcome)
    png = OUT_INDIVIDUAL / f"forest_{base}.png"
    pdf = OUT_INDIVIDUAL / f"forest_{base}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def main() -> None:
    df = load_poolable_data(DATA)
    summary = summarise_by_outcome(df)
    summary["Ratio_95CI"] = summary.apply(ratio_text, axis=1)
    summary["Prediction_interval"] = summary.apply(pi_text, axis=1)
    summary.to_csv(OUT_TABLES / "python_reml_outcome_summary.csv", index=False)

    plot_summary_forest(summary)
    for outcome in summary["Outcome"].tolist():
        plot_individual_forest(df, outcome)

    # Write a compact validation file for GitHub audit.
    audit = {
        "n_rows_poolable": len(df),
        "n_outcomes": summary["Outcome"].nunique(),
        "n_studies": df["Study_ID"].nunique() if "Study_ID" in df else df["Author_Year"].nunique(),
    }
    pd.Series(audit).to_csv(OUT_TABLES / "python_pipeline_audit.csv")
    print("Pipeline complete")
    print(audit)


if __name__ == "__main__":
    main()
