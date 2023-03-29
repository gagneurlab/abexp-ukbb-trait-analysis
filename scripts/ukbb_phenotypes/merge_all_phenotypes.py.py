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
from rep.notebook_init import setup_plot_style
setup_plot_style()

# %%
import os
import numpy as np
import pandas as pd

import json
import yaml

import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f

# import glow

# %% {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark(enable_glow=False)

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
        rule_name = 'merge_all_phenotypes',
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
meta_df = pd.read_parquet(snakemake.input["latest_meta_pq"])
meta_df

# %%
data_dfs = []

# %%
for path, df in meta_df.groupby("data_path"):
    data_df = (
        spark.read.parquet(path)
        .select([
            "eid",
            *[f.col(c) for c in df["col.name"].tolist()],
        ])
    )
    data_dfs.append(data_df)

# %%
from functools import reduce

# %%
joined_data_df = reduce(lambda df1, df2: df1.join(df2, on="eid", how="outer"), data_dfs)
joined_data_df.columns

# %%
snakemake.output

# %%
(
    joined_data_df
    .sort(["eid"])
    .write.parquet(snakemake.output["latest_data_pq"], mode="overwrite")
)

# %%
joined_data_df = spark.read.parquet(snakemake.output["latest_data_pq"])

# %%
joined_data_df.select("diagnoses_icd9_f41271_0_0").distinct().toPandas()

# %%
joined_data_df.select(meta_df.query("`field.showcase` == '41271'")["col.name"].tolist()).distinct().toPandas()

# %%
meta_df.query("`field.showcase` == '41271'")

# %%
