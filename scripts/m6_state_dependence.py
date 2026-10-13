"""M6 / A7: state dependence of the feedback parameters between 2xCO2 and 4xCO2.

(1) beta, gamma (COU-BGC and RAD) at 2x and 4x per model (esm-deck, piControl baseline), and upstream's
    land values (TCREsource_betagamma.csv) for comparison;
(2) eps scaling: eps / (dCa * dT_COU) at 2x and 4x (constant if eps ~ (d beta / dT) dCa dT);
(3) transfer test: TCRE and AF at 4x reconstructed from the 2x parameters (alpha taken at 4x), against
    the actual 4x values, under each convention.
Writes results/m6_state_dependence.csv, results/m6_state_dependence.json, results/fig_m6_state_dependence.png
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tcre_evidence import feedbacks as fb
from tcre_evidence import identity as idt
from tcre_evidence import plotting as pl
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS


def main():
    pl.style()
    t2, t4 = fb.table("2x", "picontrol"), fb.table("4x", "picontrol")
    D = pd.DataFrame(index=t2.index)
    for v in ["beta_land_coubgc", "beta_ocean_coubgc", "gamma_land_coubgc", "gamma_ocean_coubgc", "gamma_land_rad",
              "gamma_ocean_rad", "eps_land", "eps_ocean", "tcre", "af", "alpha"]:
        D[f"{v}_2x"], D[f"{v}_4x"] = t2[v], t4[v]
    for x in ("land", "ocean"):
        D[f"eps_{x}_scaled_2x"] = t2[f"eps_{x}"] / (t2.dCa * t2.dT_cou)
        D[f"eps_{x}_scaled_4x"] = t4[f"eps_{x}"] / (t4.dCa * t4.dT_cou)
    up2, up4 = ud.land_csv("2xCO2"), ud.land_csv("4xCO2")
    D["gamma_land_upcsv_2x"] = up2.gamma_land.reindex(model=D.index).values
    D["gamma_land_upcsv_4x"] = up4.gamma_land.reindex(model=D.index).values
    js = {}
    for v in ["gamma_land_coubgc", "gamma_ocean_coubgc", "gamma_land_rad", "gamma_ocean_rad", "beta_land_coubgc",
              "beta_ocean_coubgc", "gamma_land_upcsv"]:
        d = D[f"{v}_4x"] - D[f"{v}_2x"]
        js[v] = dict(mean_change=float(d.mean()), n_more_negative=int((d < 0).sum()), n=int(d.notna().sum()),
                     mean_2x=float(D[f"{v}_2x"].mean()), mean_4x=float(D[f"{v}_4x"].mean()))
    for x in ("land", "ocean"):
        js[f"eps_{x}_scaled"] = dict(mean_2x=float(D[f"eps_{x}_scaled_2x"].mean()), mean_4x=float(D[f"eps_{x}_scaled_4x"].mean()),
                                    corr_2x_4x=float(D[f"eps_{x}_scaled_2x"].corr(D[f"eps_{x}_scaled_4x"])))
    # transfer: 2x parameters, 4x alpha
    tr = {}
    for conv, (b, g) in {"coubgc": ("coubgc", "coubgc"), "bgcrad": ("bgcrad", "rad")}.items():
        beta = t2[f"beta_land_{b}"] + t2[f"beta_ocean_{b}"]
        gam = t2[f"gamma_land_{g}"] + t2[f"gamma_ocean_{g}"]
        rec = idt.tcre(t4.alpha, beta, gam)
        err = rec - t4.tcre
        D[f"tcre4x_from2x_{conv}"] = rec
        tr[conv] = dict(err_mean=float(err.mean()), err_sd=float(err.std()), max_abs=float(err.abs().max()),
                        worst=err.abs().idxmax(), rel_err_mean=float((err / t4.tcre).mean()))
    js["transfer_2x_to_4x"] = tr
    D.to_csv(RESULTS / "m6_state_dependence.csv", float_format="%.4g")
    json.dump(js, open(RESULTS / "m6_state_dependence.json", "w"), indent=1)
    print(json.dumps(js, indent=1))

    fig, axes = plt.subplots(1, 3, figsize=(9, 2.9))
    for ax, v, lab in [(axes[0], "gamma_land_coubgc", "γ$_L$ COU−BGC"), (axes[1], "gamma_ocean_coubgc", "γ$_O$ COU−BGC")]:
        x, y = D[f"{v}_2x"], D[f"{v}_4x"]
        ax.plot(x, y, "o", ms=5, color=pl.SERIES[0])
        lim = [min(x.min(), y.min()) * 1.1, max(0, x.max(), y.max()) + 5]
        ax.plot(lim, lim, color=pl.TEXT2, lw=0.8, ls=":")
        for m in D.index:
            ax.annotate(m.split("-")[0], (x[m], y[m]), fontsize=6, color=pl.TEXT2, xytext=(3, 2), textcoords="offset points")
        ax.set_xlabel(f"{lab} at 2×CO$_2$ (PgC K$^{{-1}}$)")
        ax.set_ylabel(f"at 4×CO$_2$")
    ax = axes[2]
    for j, conv in enumerate(["coubgc", "bgcrad"]):
        ax.plot(D.tcre_4x, D[f"tcre4x_from2x_{conv}"], "o", ms=5, color=pl.SERIES[j],
                label={"coubgc": "2× COU−BGC parameters", "bgcrad": "2× BGC/RAD parameters"}[conv])
    lim = [D.tcre_4x.min() * 0.9, D.tcre_4x.max() * 1.1]
    ax.plot(lim, lim, color=pl.TEXT2, lw=0.8, ls=":")
    ax.set_xlabel("actual TCRE at 4×CO$_2$ (K EgC$^{-1}$)")
    ax.set_ylabel("reconstructed from 2× β, γ")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    pl.save(fig, RESULTS / "fig_m6_state_dependence.png")


if __name__ == "__main__":
    main()
