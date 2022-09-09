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
spark = init_spark()
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

print(json.dumps(snakemake.__dict__, indent=2, default=str))

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

# # Save output

snakemake.output['data_pq']

(
    df
    .sort(["gene", "individual"])
#     .repartitionByRange(128, ["gene", "individual"])
    .write
    .parquet(snakemake.output['data_pq'], mode="overwrite")
)

spark.stop()


