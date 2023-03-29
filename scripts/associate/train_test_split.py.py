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
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'train_test_split',
        default_wildcards={
            "phenotype_col": "c_reactive_protein",
            # "phenotype_col": "systolic_blood_pressure",
            # "phenotype_col": "systolic_blood_pressure_f4080_0_0",
            "covariates": "sex_age_genPC_CLMP_PRS",
            # "covariates": "sex+age+genPC+CLMP",
            # "covariates": "sex_age_genPC",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %% [markdown]
# # Load configuration

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %% {"tags": []}
covariates_df = (
    pl.scan_parquet(snakemake.input["covariates_pq"])
    .select([
        "individual",
        phenotype_col,
    ])
)
covariates_df.schema

# %% {"tags": []}
samples_df = covariates_df.collect().to_pandas()
samples_df

# %%
is_boolean_phenotype = covariates_df.schema[phenotype_col] == pl.Boolean()
is_boolean_phenotype

# %% [markdown]
# # train/test split

# %%
import sklearn.model_selection

# %% {"tags": []}
snakemake.params

# %%
if is_boolean_phenotype:
    # train/test
    # train_test_split = sklearn.model_selection.StratifiedKFold(
    #     n_splits=2,
    #     random_state=42,
    #     shuffle=True,
    # )
    
    # data for association testing
    association_split = sklearn.model_selection.StratifiedShuffleSplit(
        train_size=snakemake.params["assocation_split"],
        random_state=42,
    )
    # split train set into validation groups
    valid_split = sklearn.model_selection.StratifiedKFold(
        n_splits=snakemake.params["phenotype_prediction_folds"],
        random_state=42,
        shuffle=True,
    )
else:
    # # train/test
    # train_test_split = sklearn.model_selection.KFold(
    #     n_splits=2,
    #     random_state=42,
    #     shuffle=True,
    # )
    
    # data for association testing
    association_split = sklearn.model_selection.ShuffleSplit(
        train_size=snakemake.params["assocation_split"],
        random_state=42,
    )
    # split train set into validation groups
    valid_split = sklearn.model_selection.KFold(
        n_splits=snakemake.params["phenotype_prediction_folds"],
        random_state=42,
        shuffle=True,
    )


# %%
samples_df["fold"] = ""
association_index, train_index = next(association_split.split(samples_df["individual"], samples_df[phenotype_col]))
samples_df.loc[association_index, "fold"] = "association_testing"

for idx, (train_fold_index, test_fold_index) in enumerate(valid_split.split(train_index, samples_df[phenotype_col].iloc[train_index])):
    print(idx)
    print(test_fold_index)
    
    fold_index = train_index[test_fold_index]
    samples_df.loc[fold_index, "fold"] = f"fold {idx + 1}"

# %%
samples_df["fold"].value_counts().sort_index()

# %%
samples_df

# %% [markdown]
# # write output

# %%
snakemake.output

# %%
samples_df.to_csv(snakemake.output["train_test_split_tsv"], sep="\t", index=False)

# %%
samples_df.to_parquet(snakemake.output["train_test_split_pq"], index=False)

# %%
