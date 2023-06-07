# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian3]
#     language: python
#     name: conda-env-anaconda-florian3-py
# ---

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
# # code autoreload
# # %load_ext autoreload
# # %autoreload 2
import os
import sys
import shutil

import json
import yaml

from pprint import pprint

import collections
import random
import math
import numpy as np
import pandas as pd

from cached_property import cached_property
# import joblib

import xarray as xr
import dask
import dask.dataframe as ddf
import dask.array as da
import zarr

try:
    from tqdm.notebook import tqdm
except:
    from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib

# %matplotlib inline
# %config InlineBackend.figure_format='retina'

import seaborn as sns
import plotnine as pn

# %%
from rep.notebook_init import init_ray, setup_plot_style
init_ray()
setup_plot_style()

# %%
import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f

import rep.spark_functions as sf

from rep.notebook_init import init_spark
spark = init_spark()
spark

# %%
import logging
logging.basicConfig()


# %%
snakefile_path = os.getcwd() + "/../../../Snakefile"

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'feature_sets__expr',
        default_wildcards={
#             'ds_dir': 'gtex_noCAT_samplefilter_maxol20_privvar',
            'ds_dir': 'gtex_noCATskin_samplefilter_maxol20_privvar',
#             'ds_dir': 'prokisch',
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %% [markdown]
# # Load input data

# %%
with open(snakemake.input["ds_config"], "r") as fd:
    config = yaml.safe_load(fd)
print(json.dumps(config, indent=2))

# %%
filter_gene_min_mu = config.get("filter_gene_min_mu", 0)
if filter_gene_min_mu is None:
    filter_gene_min_mu = 0
filter_gene_min_mu

# %%
expression_xrds = xr.open_zarr(snakemake.input["outrider_xrds"])
expression_xrds

# %%
expression_df = spark.read.parquet(snakemake.input["outrider_pq"])
expression_df.printSchema()

# %%
variables = [
    ('zscore', 0),
    ('pval', 0),
    ('padj', 0),
    ('hilo_padj', 0),
    ('l2fc', 0),
    ('counts', 0),
    ('missing', True),
    ('mu', 0),
    ('theta', 0),
]
variables

# %%
assert all([v in expression_df.columns for v, default_value in variables])

# %% [markdown]
# # Convert to features

# %%
expression_df.select("subtissue").distinct().sort("subtissue").toPandas()

# %%
# aggregated_df.select(
#         "individual",
#         "gene",
#         # rename all columns to url-encoded versions in order to escape special characters
#         *[f.col(f"`{c}`").alias(quote(c)) for c in expression_xrds.subtissue.values]
#     )

# %%
from urllib.parse import quote, unquote

aggregated_df = (
    expression_df
    .groupby("individual", "gene")
    .pivot("subtissue")#, expression_xrds.subtissue.values)
    .agg(
        f.struct(
            *[
                f.coalesce(
                    f.first(f.col(c)),
                    f.lit(default)
                ).alias(c) for c, default in variables
            ],
            (f.first(f.col("mu")) >= f.lit(filter_gene_min_mu)).alias("is_OUTRIDER_callable")
        )
    )
    .select(
        "individual",
        "gene",
        # rename all columns to url-encoded versions in order to escape special characters
        *[f.col(c).alias(sf.normalise_name(c)) for c in expression_xrds.subtissue.values]
    )
)
aggregated_df.printSchema()

# %%
snakemake.output['data_pq']

# %%
(
    aggregated_df
    .sort("individual", "gene")
    .write
    .parquet(snakemake.output['data_pq'], mode="overwrite")
)

# %%
