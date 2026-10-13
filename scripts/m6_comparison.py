"""M6: comparison of TCRE lines of evidence (HANDOVER §5.4), as interval bars (5-17-50-83-95 %).

(i)   raw CMIP6: esm-deck 1pctCO2 TCRE (11 models, piControl baseline)
(ii)  upstream framework as published (M3 rerun, variant U)
(iii) upstream with joint land model and TCR-carbon rank correlation from the deck (variant J+a)
(iv)  TCR x R, R unconstrained (TCR ~ N(1.7, 0.36), as upstream)
(v)   TCR x R, R constrained by the 1960-2014 ocean sink (same TCR prior)
(vi)  historical attribution (upstream TCRE_from_DAMIP.ipynb saved output: 1.386 +/- 0.087; quoted, not rerun)
Perfect-model CRPS skill against the raw ensemble is annotated where it exists (A4 for ii; M5 for v).
Writes results/fig_m6_comparison.png, results/m6_comparison.csv
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tcre_evidence import feedbacks as fb
from tcre_evidence import plotting as pl
from tcre_evidence.paths import RESULTS

Q = [5, 17, 50, 83, 95]


def main():
    pl.style()
    rng = np.random.default_rng(0)
    m3 = np.load(RESULTS / "m3_samples.npz")
    m5 = np.load(RESULTS / "m5_samples.npz")
    k5 = {k: k for k in m5.files}
    key = lambda mode, h: [k for k in m5.files if k.startswith("upstream_N1.7_0.36__R__") and k.endswith(f"__{mode}__{h}")][0]
    lines = {
        "(i) raw CMIP6 1pctCO2 (11 models)": fb.table("2x", "picontrol").tcre.values,
        "(ii) process framework as published": m3["U"],
        "(iii) framework, joint land + TCR–carbon covariance": m3["J_a"],
        "(iv) TCR × R, R unconstrained": m5[key("unconstrained", "-")],
        "(v) TCR × R, R | ocean sink 1960–2014": m5[key("substitute", "H_ocean")],
        "(vi) historical attribution (DAMIP, quoted)": rng.normal(1.386, 0.087, 20000),
    }
    skill = {}
    try:
        s4 = json.load(open(RESULTS / "m4_summary.json"))["skill"]
        sk = [s for s in s4 if s["mode"] == "loo" and s["quantity"] == "tcre"][0]
        skill["(ii) process framework as published"] = f"perfect-model CRPSS {sk['crpss_vs_raw']:+.2f} ({sk['n_better']}/{sk['n']} better)"
    except Exception:
        pass
    a5 = json.load(open(RESULTS / "m5.json"))["perfect_model"]
    v = [r for r in a5 if r["target"] == "R" and r["H"] == "H_ocean" and not r["lfo"] and r["obs_err"] == 1.0
         and r["mode"] == "H_only"][0]
    skill["(v) TCR × R, R | ocean sink 1960–2014"] = f"perfect-model CRPSS {v['crpss_vs_ensemble']:+.2f}"
    rows = []
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    for i, (lab, s) in enumerate(lines.items()):
        q = np.percentile(s, Q)
        rows.append(dict(line=lab, mean=float(np.mean(s)), **{f"p{p:02d}": float(x) for p, x in zip(Q, q)}))
        col = pl.SERIES[0] if i else pl.TEXT2
        ax.hlines(i, q[0], q[4], color=col, lw=1.5)
        ax.hlines(i, q[1], q[3], color=col, lw=5)
        ax.plot(q[2], i, "o", ms=6, color=pl.SURFACE, mec=col, mew=2)
        if i == 0:
            ax.plot(s, np.full_like(s, i) + 0.18, "|", ms=6, color=pl.TEXT2)
        if lab in skill:
            ax.annotate(skill[lab], (q[4], i), xytext=(6, -3), textcoords="offset points", fontsize=7, color=pl.TEXT2)
    ax.set_yticks(range(len(lines)), list(lines))
    ax.invert_yaxis()
    ax.set_xlabel("TCRE (K EgC$^{-1}$); bars 5–95 % (thin) and 17–83 % (thick), circle median")
    ax.set_xlim(0.6, 3.6)
    pl.save(fig, RESULTS / "fig_m6_comparison.png")
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "m6_comparison.csv", index=False, float_format="%.3f")
    print(df.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
