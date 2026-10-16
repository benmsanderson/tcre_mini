"""M10: does the development configuration's own nitrogen decomposition reproduce the end-to-end difference?

Fits the configuration's beta_L and gamma_L process models (default priors; gamma_L also with the paired-
experiment nu evidence as in m6_identifiability) and, for each posterior draw, predicts each esm-deck model's
expected beta_L and gamma_L from its own process flags (baseline x prod eta^L; baseline + sum delta L, tropics
+ extratropics; no model-level random effects), then its carbon factor R at 2xCO2 (dCa 285 ppm, TCR 1.8 K,
ocean at the configuration's posterior means). Compares
  (1) the implied ratio of group-mean R, nitrogen-flagged over unflagged models, with the observed ratio from
      the deck 1pctCO2 runs (results/m7_nitrogen_historical.json, 1.063);
  (2) the implied pure nitrogen effect on R (switch the flag on for every model, other flags as they are);
  (3) the same with the gamma_L nitrogen term removed (beta route only).
Writes results/m10_nitrogen_decomposition.json; with --beta-split, results/m10_beta_split.json (the
'all processes on' beta_L split into baseline, nitrogen factor and other factors)
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import json

import numpy as np
import pandas as pd
import pymc as pm
from joblib import Parallel, delayed

from tcre_evidence import feedbacks as fb
from tcre_evidence import upstream as up
from tcre_evidence import upstream_data as ud
from tcre_evidence.paths import RESULTS, SNAPSHOT_DATA
from tcre_evidence.upstream_patched import model_builders as mb

SEED = 0
SAMPLER = dict(draws=1000, tune=1000, chains=4, cores=1, target_accept=0.99, progressbar=False)
DCA, TCR, K = 285.0, 1.8, 2.12


def nu_evidence():
    g = pd.read_csv(SNAPSHOT_DATA / "TCREsource_betagamma.csv").set_index("model").gamma_L_2xCO2
    return {"nitrogen": [float(g["ACCESS-ESM_CN"] - g["ACCESS-ESM_C"]), float(g["UKESM1-1_ctrl"] - g["UKESM1-1_nonlim"])],
            "fire": [float(g["UKESM1-1_fire"] - g["UKESM1-1_ctrl"])], "veg": [float(g["UKESM1-1_ctrl"] - g["UKESM1-1_nodgvm"])]}


def post(idata, v):
    x = idata.posterior[v].values
    return x.reshape(x.shape[0] * x.shape[1], *x.shape[2:])


def fit(kind):
    inp = up.inputs()
    if kind == "beta":
        d = inp["beta_land_data"]
        m = mb.build_beta_land_model(d["beta_land"], d["lookup_table"], up.priors.beta_land)
        with m:
            idata = pm.sample(random_seed=SEED + 1, **SAMPLER)
        return kind, dict(U=post(idata, "beta_land_unscaled"), log_eta=post(idata, "beta_land_log_eta"),
                          procs=list(up.priors.beta_land["eta"].keys())), up.diagnostics(idata, ["beta_land"])
    d = inp["gamma_land_data"]
    m = mb.build_gamma_land_model(d["gamma_land_tropics"], d["gamma_land_extratropics"], d["lookup_table"],
                                  inp["iav_df"], inp["iav_obs_mean"], inp["iav_obs_sigma"], up.priors.gamma_land)
    procs = list(d["lookup_table"].process.values)
    if kind == "gamma_nu":
        with m:
            for p, vals in nu_evidence().items():
                i = procs.index(p)
                pm.Normal(f"nu_{p}", mu=m["gamma_land_tropics_delta"][i] + m["gamma_land_extratropics_delta"][i],
                          sigma=5.0, observed=np.array(vals))
    with m:
        idata = pm.sample(random_seed=SEED + 2, **SAMPLER)
    U = post(idata, "gamma_land_tropics_unscaled") + post(idata, "gamma_land_extratropics_unscaled")
    D = post(idata, "gamma_land_tropics_delta") + post(idata, "gamma_land_extratropics_delta")
    return kind, dict(U=U, delta=D, procs=procs), up.diagnostics(idata, ["gamma_land"])


def R_of(beta_L, gamma_L, bO, gO):
    return 1e3 / (K * DCA + (beta_L + bO) * DCA + (gamma_L + gO) * TCR)


def main():
    res = Parallel(n_jobs=3)(delayed(fit)(k) for k in ("beta", "gamma", "gamma_nu"))
    F = {k: v for k, v, _ in res}
    diag = {k: d for k, _, d in res}
    oc = up.fit_ocean(up.inputs(), seed=SEED + 3, sampler=dict(cores=1))
    bO, gO = float(up.flat(oc.idata, "beta_ocean").mean()), float(up.flat(oc.idata, "gamma_ocean").mean())
    diag["ocean"] = oc.diagnostics
    t = fb.table("2x", "picontrol")
    models = [m for m in t.index]
    lt = ud.lookup_table().sel(model=models)
    Lb = np.array([[float(lt.sel(model=m, process=p)) for p in F["beta"]["procs"]] for m in models])
    obs = json.load(open(RESULTS / "m7_nitrogen_historical.json"))["R"]
    obs_ratio = obs["mean_with_N"] / obs["mean_without_N"]
    out = dict(observed_group_ratio=obs_ratio, observed_perm_p=obs["perm_p"], ocean=dict(beta=bO, gamma=gO),
               diagnostics=diag, settings=dict(dCa=DCA, TCR=TCR))
    nb = F["beta"]["procs"].index("nitrogen")
    for gk in ("gamma", "gamma_nu"):
        Lg = np.array([[float(lt.sel(model=m, process=p)) for p in F[gk]["procs"]] for m in models])
        ng = F[gk]["procs"].index("nitrogen")
        n = min(len(F["beta"]["U"]), len(F[gk]["U"]))
        U_b, le = F["beta"]["U"][:n], F["beta"]["log_eta"][:n]
        U_g, dl = F[gk]["U"][:n], F[gk]["delta"][:n]

        def predict(Lbm, Lgm, gamma_N=True):
            beta = U_b[:, None] * np.exp(le @ Lbm.T)
            dg = dl.copy()
            if not gamma_N:
                dg[:, ng] = 0.0
            gamma = U_g[:, None] + dg @ Lgm.T
            return R_of(beta, gamma, bO, gO)  # draws x models

        isN = Lb[:, nb] == 1
        R = predict(Lb, Lg)
        grp = R[:, isN].mean(1) / R[:, ~isN].mean(1)
        Lb1, Lg1, Lb0, Lg0 = Lb.copy(), Lg.copy(), Lb.copy(), Lg.copy()
        Lb1[:, nb] = 1; Lg1[:, ng] = 1; Lb0[:, nb] = 0; Lg0[:, ng] = 0
        pure = (predict(Lb1, Lg1) / predict(Lb0, Lg0)).mean(1)
        pure_beta_only = (predict(Lb1, Lg1, gamma_N=False) / predict(Lb0, Lg0, gamma_N=False)).mean(1)
        q = lambda x: dict(mean=float(x.mean()), p05=float(np.percentile(x, 5)), p95=float(np.percentile(x, 95)),
                           prob_gt_observed=float((x > obs_ratio).mean()))
        out[gk] = dict(group_ratio=q(grp), pure_N_effect=q(pure), pure_N_effect_beta_only=q(pure_beta_only),
                       eta_N=float(np.exp(le[:, nb]).mean()), delta_N_total=float(dl[:, ng].mean()),
                       delta_N_total_sd=float(dl[:, ng].std()))
    json.dump(out, open(RESULTS / "m10_nitrogen_decomposition.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != "diagnostics"}, indent=1, default=float))
    print({k: (round(v["rhat_max"], 4), v["divergences"]) for k, v in diag.items()})




def beta_split():
    """Split the configuration's 'all processes on' beta_L into baseline, nitrogen factor and other factors."""
    _, f, d = fit("beta")
    eta = np.exp(f["log_eta"])
    procs = f["procs"]
    U = f["U"]
    sh = ud.cmip_sheet_dataset().drop_dims("process").to_dataframe()
    lt = ud.lookup_table().to_pandas().T
    b = sh.beta_land.dropna()
    nflag = lt.reindex(b.index).nitrogen == 1
    out = dict(U=float(U.mean()), eta={p: float(eta[:, i].mean()) for i, p in enumerate(procs)},
               beta_allon=float((U * eta.prod(1)).mean()),
               beta_U_times_etaN=float((U * eta[:, procs.index("nitrogen")]).mean()),
               cmip_mean=float(b.mean()), cmip_mean_N=float(b[nflag].mean()), cmip_mean_noN=float(b[~nflag].mean()),
               diagnostics=d)
    json.dump(out, open(RESULTS / "m10_beta_split.json", "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    import sys
    beta_split() if "--beta-split" in sys.argv else main()
