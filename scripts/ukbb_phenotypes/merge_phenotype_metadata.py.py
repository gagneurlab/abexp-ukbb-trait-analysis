# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

# %%
# # code autoreload
# # %load_ext autoreload
# # %autoreload 2
import os
import sys
import shutil

from pprint import pprint

import collections
import random
import math
import numpy as np
import pandas as pd

import json
import yaml

# import joblib

try:
    from tqdm.notebook import tqdm
except:
    from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib

# import seaborn as sns
# import plotnine as pn




# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

# %%
# from rep.notebook_init import init_ray, setup_plot_style
# init_ray()
# setup_plot_style()

# %%
snakefile_path = os.getcwd() + "/../Snakefile"
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
        rule_name = 'merge_phenotype_metadata',
        default_wildcards={
#             'ds_dir': 'noCAT_samplefilter_maxol20_privvar'
#             'ds_dir': 'gtex_noCAT_samplefilter_maxol20_privvar'
#             'ds_dir': 'gtex_noCATskin_samplefilter_maxol20_privvar'
#             'ds_dir': 'gtex_v8_old_dna'
#             'ds_dir': 'prokisch_WGS'
#             'ds_dir': 'prokisch_WGS_maxol150'
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %%
meta_paths_df = pd.DataFrame(snakemake.params["path_codes"], columns=["meta_pq", "dirname", "ukbb_code"]).sort_values("ukbb_code", ascending=False)
meta_paths_df

# %%
meta_dfs = []

for idx, row in meta_paths_df.iterrows():
    meta_sub_df = pd.read_parquet(row["meta_pq"])
    # drop 'eid' description
    meta_sub_df = meta_sub_df.iloc[1:]
    meta_sub_df["data_path"] = os.path.join(os.path.dirname(row["meta_pq"]), f"{row['ukbb_code']}.data.parquet")
    meta_sub_df["ukbb_code"] = row["ukbb_code"]
    
    meta_dfs.append(meta_sub_df)
meta_df = pd.concat(meta_dfs, axis=0)

# %%
meta_df

# %%
is_latest_df = (
    meta_df[["field.tab", "ukbb_code"]]
    .sort_values("ukbb_code", ascending=False)
    .groupby("field.tab")
    .first()
)
is_latest_df["latest"] = True
is_latest_df

# %%
meta_df_annotated = (
    meta_df
    .merge(is_latest_df, on=["field.tab", "ukbb_code"], how="left")
    .fillna({"latest": False})
)
meta_df_annotated

# %%
meta_df_annotated.to_parquet(snakemake.output["merged_meta_pq"])

# %%
meta_df_annotated.query("`latest` == True").to_parquet(snakemake.output["latest_meta_pq"])

# %%
