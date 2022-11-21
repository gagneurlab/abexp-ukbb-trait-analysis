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

import glow

# %% {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark()

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
        rule_name = 'decode_phenotype',
        default_wildcards={
            "phenotype": "severe_LDL",
            # "phenotype": "Diabetes",
            # "phenotype": "standing_height",
            # "phenotype": "Asthma",
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))

# %%
meta_df = pd.read_parquet(snakemake.input["latest_meta_pq"])
meta_df

# %%
with open(snakemake.input["phenotype_coding_yaml"], "r") as fd:
    phenotype_coding = yaml.safe_load(fd)
print(json.dumps(phenotype_coding, indent=2, default=str))

# %% [markdown]
# # find required fields

# %% {"tags": []}
# schema:
# { <aggregation>: {<col_regex>: [<list of matching terms>]} }
required_columns = set()
for regex in phenotype_coding["fields"]:
    matching_cols = meta_df["col.name"].str.match(regex)
    matching_colnames = meta_df["col.name"].loc[matching_cols].tolist()
    required_columns.update(matching_colnames)
required_columns

# %% [markdown] {"tags": []}
# # read data

# %%
data_dfs = []

# %%
for path, df in meta_df[meta_df["col.name"].isin(required_columns)].groupby("data_path"):
    data_df = (
        spark.read.parquet(path)
        .select([
            "eid",
            *[f.col(c) for c in df["col.name"].tolist()],
        ])
    )
    data_dfs.append(data_df)

# %%
len(data_dfs)

# %%
from functools import reduce

# %% {"tags": []}
joined_data_df = reduce(lambda df1, df2: df1.join(df2, on="eid", how="outer"), data_dfs)
joined_data_df.columns

# %% [markdown]
# # aggregate columns

# %%
if phenotype_coding["type"] == "any":
    # schema:
    # { <aggregation>: {<col_regex>: [<list of matching terms>]} }
    # expand terms
    col_terms = []
    for regex, terms in phenotype_coding["fields"].items():
        matching_cols = meta_df["col.name"].str.match(regex)
        matching_colnames = meta_df["col.name"].loc[matching_cols]

        for col in matching_colnames:
            # expr = f.coalesce(f.col(col).isin(terms), f.lit(False))
            expr = f.col(col).isin(terms)
            col_terms.append(expr)
    col_terms

    phenotype_expr = reduce(lambda a, b: a | b, col_terms)


# %%
def mean_agg(col, ignore_null=True, ignore_nan=True):
    if ignore_null:
        col = f.filter(col, lambda x: ~ f.isnull(x))
    if ignore_nan:
        col = f.filter(col, lambda x: ~ f.isnan(x))
    
    def merge(acc, x):
        count = acc.count + 1
        sum = acc.sum + x
            
        return f.struct(
            count.alias("count"),
            sum.alias("sum")
        )

    return f.aggregate(
        col,
        initialValue=f.struct(
            f.lit(0).alias("count"), 
            f.lit(0.0).alias("sum"),
        ),
        merge=merge,
        finish=lambda acc: acc.sum / acc.count,
    ).alias("mean")


# %%
def sum_agg(col, ignore_null=True, ignore_nan=True):
    if ignore_null:
        col = f.filter(col, lambda x: ~ f.isnull(x))
    if ignore_nan:
        col = f.filter(col, lambda x: ~ f.isnan(x))
    
    def merge(acc, x):
        sum = acc + x
        return sum

    return f.aggregate(
        col,
        initialValue=f.lit(0.0),
        merge=merge,
    ).alias("sum")


# %%
def median_agg(col, ignore_null=True, ignore_nan=True):
    if ignore_null:
        col = f.filter(col, lambda x: ~ f.isnull(x))
    if ignore_nan:
        col = f.filter(col, lambda x: ~ f.isnan(x))
    
    # sort array
    col = f.array_sort(col)
    
    return f.when(
        f.size(col) % 2 == 0,
        # mean of two middle values
        (col[f.floor(f.size(col) / 2)] + col[f.ceil(f.size(col) / 2)]) / 2
    ).otherwise(
        col[f.floor(f.size(col) / 2)]
    ).alias("median")


# %%
if phenotype_coding["type"] == "mean":
    phenotype_expr = mean_agg(f.array(matching_colnames))

if phenotype_coding["type"] == "median":
    phenotype_expr = median_agg(f.array(matching_colnames))

if phenotype_coding["type"] == "sum":
    phenotype_expr = sum_agg(f.array(matching_colnames))

# %%
# joined_data_df.select("eid", median_agg(f.array(matching_colnames)), mean_agg(f.array(matching_colnames)), *matching_colnames).toPandas().dropna()

# %% [markdown]
# ## create dataframe

# %%
phenotype_expr = phenotype_expr.alias(snakemake.wildcards["phenotype"])
phenotype_expr

# %%
phenotype_df = joined_data_df.select([
    "eid",
    phenotype_expr,
])
phenotype_df.printSchema()

# %% [markdown]
# ## apply expression if existing

# %%
phenotype_coding["expression"].format(col=snakemake.wildcards["phenotype"])

# %%
f.expr(phenotype_coding["expression"].format(col=snakemake.wildcards["phenotype"])).alias(snakemake.wildcards["phenotype"])

# %%
phenotype_df = phenotype_df.withColumn(snakemake.wildcards["phenotype"],
    f.expr(phenotype_coding["expression"].format(col="`" + snakemake.wildcards["phenotype"] + "`"))
)
phenotype_df.printSchema()

# %%
snakemake.output

# %%
(
    phenotype_df
    .sort(["eid"])
    .write.parquet(snakemake.output["phenotype_pq"], mode="overwrite")
)

# %%
phenotype_df = spark.read.parquet(snakemake.output["phenotype_pq"])

# %%
phenotype_pd_df = phenotype_df.toPandas()
phenotype_pd_df

# %%
phenotype_df.

# %%
pyspark.sql.types.to

# %%

# %%
pyspark.sql.types.

# %%
x.dataType.simpleString()

# %%

# %%
import plotnine as pn

# %%
pheno_dtype = phenotype_df.schema[snakemake.wildcards["phenotype"]].dataType
pheno_dtype

# %%
if pheno_dtype == t.BooleanType():
    plot = (
        pn.ggplot(phenotype_pd_df, pn.aes(x=snakemake.wildcards["phenotype"]))
        + pn.geom_bar(stat="count")
    )
elif (
    pheno_dtype == t.ShortType(),
    pheno_dtype == t.IntegerType(),
    pheno_dtype == t.LongType(),
    pheno_dtype == t.FloatType(),
    pheno_dtype == t.DoubleType(),
):
    plot = (
        pn.ggplot(phenotype_pd_df, pn.aes(x=snakemake.wildcards["phenotype"]))
        + pn.geom_histogram()
    )
else:
    plot = None

# %%
plot

# %%
if plot is not None:
    plot.save(snakemake.params["output_basepath"] + "/data_plot.pdf", dpi=450)
    plot.save(snakemake.params["output_basepath"] + "/data_plot.png", dpi=450)

# %%
