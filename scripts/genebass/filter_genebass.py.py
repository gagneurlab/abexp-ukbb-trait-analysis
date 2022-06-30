# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

from IPython.display import display

# +
import os
import numpy as np
import pandas as pd

import json
import yaml

import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f

import glow


# + {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark()
# -

spark

snakefile_path = os.getcwd() + "/../Snakefile"
snakefile_path

# +
# del snakemake
# -

try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'filter_genebass',
        default_wildcards={
        }
    )

print(json.dumps(snakemake.__dict__, indent=2, default=str))

# +
# RESULTS_MT = "/s/project/rep/processed/genebass/results.mt"
# RESULTS_TSV = "/s/project/rep/processed/genebass/results.tsv"
# RESULTS_PQ = "/s/project/rep/processed/genebass/results.parquet"

# FILTERED_PQ = "/s/project/rep/processed/genebass/filtered.parquet"

# BONFERRONI_CUTOFF = 0.01
# -

RESULTS_PQ = snakemake.input["results_pq"]
FILTERED_PQ = snakemake.output["results_filtered_pq"]
BONFERRONI_CUTOFF = snakemake.params["bonferroni_cutoff"]

genebass_df = spark.read.parquet(RESULTS_PQ)
genebass_df.printSchema()

# +
# with pd.option_context('display.max_columns', None):
#     display(
#         genebass_df.drop(
#             *[f"Nmarker_MACCate_{i}" for i in range(1, 9)],
#         ).limit(10).toPandas()
#     )
# -

with pd.option_context('display.max_columns', None):
    display(
        genebass_df.drop(
            *[f"Nmarker_MACCate_{i}" for i in range(1, 9)],
        )
        .filter(f.col("Pvalue") < 0.05)
        .filter(f.col("annotation") == f.lit("pLoF"))
#         .limit(10)
#         .toPandas()
        .count()
    )

cutoffs_df = (
    genebass_df
    .groupby("annotation")
    .count()
    .withColumn("bonferroni_cutoff", BONFERRONI_CUTOFF / f.col("count"))
    .persist()
)

cutoffs_df.toPandas()

filtered_genebass_df = (
    genebass_df
    .join(cutoffs_df, on="annotation", how="left")
    .filter(f.col("Pvalue_Burden") < f.col("bonferroni_cutoff"))
    .drop(
        *[f"Nmarker_MACCate_{i}" for i in range(1, 9)],
    )
    .persist()
)
filtered_genebass_df.printSchema()

filtered_genebass_df.select("phenocode").distinct().toPandas().sort_values("phenocode")

filtered_genebass_df.select("annotation").distinct().toPandas()

# +
# filtered_genebass_pd_df = filtered_genebass_df.toPandas()
# filtered_genebass_pd_df
# -

(
    filtered_genebass_df
    .sort([
        "annotation",
        "gene_id",
    ])
    .write
    .parquet(FILTERED_PQ, mode="overwrite")
)


