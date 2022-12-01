# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

# %%
from IPython.display import display

# %%
import os
import numpy as np
import pandas as pd

import json
import yaml

import polars as pl

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'samples',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            # "phenotype_col": "systolic_blood_pressure_f4080_0_0",
            "covariates": "sex_age_genPC_CLMP_PRS",
            # "covariates": "sex+age+genPC+CLMP",
            # "covariates": "sex_age_genPC",
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))

# %% [markdown]
# # Load configuration

# %%
with open(snakemake.input["covariates_config"], "r") as fd:
    config = yaml.safe_load(fd)

# %%
print(json.dumps(config, indent=2, default=str))

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %%
phenocode = config["phenocode"]
phenocode

# %%
restricted_formula = config["restricted_formula"]
print(restricted_formula)

# %%
pvalue_cutoff = 0.05

# %%
samples = pd.read_csv(snakemake.input["samples_txt"], header=None, names=["individual"], dtype={"individual": "string"})
samples

# %%
import sklearn.model_selection

# %%
# train/test
train_test_split = sklearn.model_selection.KFold(
    n_splits=2,
    random_state=42,
    shuffle=True,
)
valid_split = sklearn.model_selection.KFold(
    n_splits=5,
    random_state=42,
    shuffle=True,
)

samples["fold"] = ""
train_index, test_index = next(train_test_split.split(samples))
samples.loc[test_index, "fold"] = "test"

for idx, (train_fold_index, test_fold_index) in enumerate(valid_split.split(train_index)):
    print(idx)
    print(test_fold_index)
    
    fold_index = train_index[test_fold_index]
    samples.loc[fold_index, "fold"] = f"fold {idx}"

# %%
samples["fold"].value_counts()

# %%
samples

# %% [markdown]
# # write output

# %%
snakemake.output

# %%
samples.to_csv(snakemake.output["samples_tsv"], sep="\t", index=False)

# %%
samples.to_parquet(snakemake.output["samples_pq"], index=False)

# %%
