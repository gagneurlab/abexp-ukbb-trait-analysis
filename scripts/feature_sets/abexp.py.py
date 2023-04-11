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
# del snakemake
# -

try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'feature_sets__abexp',
        default_wildcards={
            'agg': "max",
            # 'agg': "mean",
            # 'agg': "median",
        }
    )

# +
# snakemake.reload()
# -

from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

os.getcwd()

# # Load input data

# ## configuration

with open(snakemake.input["agg_config"], "r") as fd:
    agg_config = yaml.safe_load(fd)
agg_config

agg = agg_config["aggregation"]
agg

exclude = agg_config.get("exclude", {})
exclude

# ## data

df = spark.read.parquet(snakemake.input["abexp_predictions"])
df.printSchema()

# exclude values
for col, values in exclude.items():
    df = df.filter(~ f.col("subtissue").isin(values))

# # aggregate over tissues

# +
agg_df = (
    df
    .filter(f.col("gene_is_expressed"))
    .groupby("gene", "individual")
)

if agg == "max":
    agg_df = agg_df.agg(
        f.max(f.col("y_pred_proba")).alias("max_AbExp")
    )
elif agg == "mean":
    agg_df = agg_df.agg(
        f.mean(f.col("y_pred_proba")).alias("mean_AbExp")
    )
elif agg == "median":
    agg_df = agg_df.agg(
        f.expr('percentile(y_pred_proba, array(0.5))')[0].alias("median_AbExp")
        # f.median(f.col("y_pred_proba")).alias("median_AbExp")
    )
else:
    raise ValueError(f"""Unknown aggregation type: '{snakemake.wildcards["agg"]}'""")

agg_df.printSchema()
# -

# # Save output

snakemake.output['data_pq']

(
    agg_df
    .sort(["gene", "individual"])
#     .repartitionByRange(128, ["gene", "individual"])
    .write
    .parquet(snakemake.output['data_pq'], mode="overwrite")
)

spark.stop()


