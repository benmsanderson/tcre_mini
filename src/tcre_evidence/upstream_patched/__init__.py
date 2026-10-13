"""Minimally patched copy of the upstream TCRE_Assessment model builders.

Source: external/TCRE_Assessment (netzeroasap/TCRE_Assessment) at commit
98c74cb8d4ccea8327cf8d325c9fe6266a2ff122, files bayes/CMIP_process.py and
bayes/model_builders.py.

At that commit the land builders do not construct:
- build_beta_land_model passes var_name= to CMIP_process.add_process_information,
  which has no such argument (TypeError);
- build_gamma_land_model nests two build_vectorized_process_model calls in one
  pm.Model, and both register unprefixed names such as "sigma_struct"
  (ValueError: variable name already exists);
- add_emergent_constraint reads model[f"{latent}_CMIP"], which is never created
  (the CMIP likelihood is named "X_obs").

The sampler output saved in notebooks/TCRE_total.ipynb lists variables such as
"beta_land_sigma_struct" and "gamma_land_tropics_rho", so the notebook was run
with a version that prefixes internal names with var_name. The patches below
restore that behaviour and change nothing else (priors, likelihoods and model
structure are untouched):

1. Internal variable names in build_vectorized_process_model are prefixed with
   f"{var_name}_" (including process deltas, e.g. "gamma_land_tropics_nitrogen").
2. The CMIP likelihood is named f"{var_name}_CMIP".
3. add_process_information accepts var_name and reads the prefixed eta/delta.
4. model_builders imports CMIP_process from this package.

`upstream.diff` in this directory is the exact diff against the pinned files.
"""
