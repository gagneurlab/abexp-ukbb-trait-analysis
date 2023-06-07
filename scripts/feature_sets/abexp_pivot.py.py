# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
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

# -


from rep.notebook_init import init_spark
spark = init_spark(enable_glow=False)

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
        rule_name = 'feature_sets__abexp_pivot',
        default_wildcards={
            # 'agg': "max",
            # 'agg': "mean",
            # 'agg': "median",
        }
    )

from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

os.getcwd()

# # Load input data

df = (
    spark.read.parquet(snakemake.input["abexp_predictions"])
    .filter(f.col("gene_is_expressed"))
)
df.printSchema()

indivs = df.select("individual").distinct().sort("individual").toPandas()["individual"]
len(indivs)

# # pivot tissues

subtissues = df.select("subtissue").distinct().sort("subtissue").toPandas()["subtissue"]
subtissues

# +
# aggregated_df.select(
#         "individual",
#         "gene",
#         # rename all columns to url-encoded versions in order to escape special characters
#         *[f.col(f"`{c}`").alias(quote(c)) for c in expression_xrds.subtissue.values]
#     )
# -

import rep.spark_functions as sf

# +
from urllib.parse import quote, unquote

aggregated_df = (
    df
    .groupby("gene", "individual")
    .pivot("subtissue", values=subtissues.to_list())
    .agg(
        f.first(f.col("y_pred_proba")).alias("AbExp")
        # f.struct(
        #     *[
        #         f.coalesce(
        #             f.first(f.col(c)),
        #             f.lit(default)
        #         ).alias(c) for c, default in variablesf
        #     ],
        #     (f.first(f.col("mu")) >= f.lit(filter_gene_min_mu)).alias("is_OUTRIDER_callable")
        # )
    )
    .select(
        "individual",
        "gene",
        # rename all columns to url-encoded versions in order to escape special characters
        f.struct(
            *[f.col(c).alias(sf.normalise_name(c)) for c in subtissues.values]
        ).alias("AbExp")
    )
)
aggregated_df.printSchema()
# -

# # Save output

snakemake.output['data_pq']

(
    aggregated_df
    .sort(["gene", "individual"])
#     .repartitionByRange(128, ["gene", "individual"])
    .write
    .parquet(snakemake.output['data_pq'], mode="overwrite")
)

spark.stop()


