# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

# +
import os
import pandas as pd

import json
import yaml

import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f



# + {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark(enable_glow=False)
# -

spark

# +
# from rep.notebook_init import init_ray
# init_ray()
# -

snakefile_path = os.getcwd() + "/../../Snakefile"

# +
#del snakemake
# -

try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'feature_sets__loftee_plof',
        default_wildcards={
        }
    )

from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

os.getcwd()

# # Load input data

df = (
    spark.read.parquet(snakemake.input["vep_predictions"])
    .select(
        "individual",
        "gene",
        f.col("features")["LoF_HC.sum"].alias("pLoF")
    )
)
df.printSchema()

plof_stats = (
    df
    .groupby("individual")
    .agg(
        f.sum(f.col("pLoF")).alias("pLoF"),
    )
    .groupby()
    .agg(
        f.sum(f.col("pLoF")).alias("pLoF_total_calls"),
        f.mean(f.col("pLoF")).alias("pLoF_average_per_individual"),
    )
    .toPandas()
)
plof_stats

variants_df = (
    spark.read.parquet(snakemake.input["vep_variants"])
    .select(
        "chrom",
        "start",
        "end",
        "ref",
        "alt",
        "gene",
        (f.col("features")["LoF_HC.sum"] > 0).alias("is_pLoF_variant")
    )
)
variants_df.printSchema()

total_num_pLoF_variants = (
    variants_df
    .groupby("chrom", "start", "end", "ref", "alt")
    .agg(
        f.max(f.col("is_pLoF_variant")).alias("is_pLoF_variant")
    )
    .filter("is_pLoF_variant")
    .count()
)
    # .select("pLoF").distinct().toPandas()

plof_stats["total_num_pLoF_variants"] = total_num_pLoF_variants
plof_stats

# # Save output

snakemake.output['stats_tsv']

plof_stats.to_csv(snakemake.output['stats_tsv'], index=False, sep="\t")

snakemake.output['data_pq']

(
    df
    .sort(["gene", "individual"])
#     .repartitionByRange(128, ["gene", "individual"])
    .write
    .parquet(snakemake.output['data_pq'], mode="overwrite")
)

spark.stop()


