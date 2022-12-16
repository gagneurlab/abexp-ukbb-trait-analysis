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
            "phenotype_col": "Asthma",
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

# %% [markdown]
# # filter samples

# %% [markdown]
# ## load annotations

# %%
phenotype_metadata_df = pl.read_parquet(snakemake.input["phenotype_metadata_pq"])
phenotype_metadata_df

# %%
data_dfs = []
for group_df in phenotype_metadata_df.filter(pl.col("field.showcase").is_in(["22006", "22011"])).groupby("data_path"):
    data_df = (
        pl.scan_parquet(group_df["data_path"].unique()[0])
        .select([
            "eid",
            *group_df["col.name"].to_list()
        ])
    )
    data_dfs.append(data_df)

# %%
len(data_dfs)

# %%
import functools
data_df = functools.reduce(lambda df1, df2: df1.join(df2, on="eid", how="outer"), data_dfs)

# %%
# rename eid col
data_df = data_df.rename({"eid": "individual"})

# %%
# make sure that sample_id is a string
data_df = data_df.with_column(pl.col("individual").cast(pl.Utf8).alias("individual"))

# %%
# join on samples
data_df = pl.from_pandas(samples).lazy().join(data_df, on="individual", how="left")

# %%
# collect relationship id's

# %%
genetic_relatedness_pairing_cols = [c for c in data_df.columns if c.startswith("genetic_relatedness_pairing_f22011")]
genetic_relatedness_pairing_cols

# %%
data_df = (
    data_df
    .with_column(
        pl.concat_list(genetic_relatedness_pairing_cols)
        .arr.eval(pl.element().drop_nulls())
        .alias("genetic_relatedness_pairing")
    )
    .drop(genetic_relatedness_pairing_cols)
)

# %%
data_df = data_df.collect()

# %%
data_df

# %% [markdown]
# ## filter for unrelated caucasians

# %%
filtered_data_df = (
    data_df
    .filter(pl.col("genetic_relatedness_pairing").arr.lengths() == 0)
    .filter(pl.col("genetic_ethnic_grouping_f22006_0_0") == pl.lit("Caucasian"))
)
filtered_data_df

# %% {"tags": []}
samples_df = filtered_data_df.select("individual").to_pandas()
samples_df

# %% [markdown]
# # write output

# %%
snakemake.output

# %%
filtered_data_df.select("individual").write_csv(snakemake.output["samples_tsv"], sep="\t")

# %%
filtered_data_df.select("individual").write_parquet(snakemake.output["samples_pq"], use_pyarrow=True, compression="snappy")

# %%
