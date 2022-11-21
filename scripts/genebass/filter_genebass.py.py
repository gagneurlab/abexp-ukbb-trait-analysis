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

# import glow


# + {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark(enable_glow=False)
# -

spark

snakefile_path = os.getcwd() + "/../../Snakefile"
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

filtered_genebass_df.select("annotation").distinct().toPandas()

# +
# filtered_genebass_pd_df = filtered_genebass_df.toPandas()
# filtered_genebass_pd_df
# -

# # write filtered genebass df

(
    filtered_genebass_df
    .sort([
        "annotation",
        "gene_id",
    ])
    .write
    .parquet(FILTERED_PQ, mode="overwrite")
)

# # Filter phenotypes

phenotype_metadata_df = spark.read.parquet(snakemake.input["phenotype_metadata_pq"])
phenotype_metadata_df = phenotype_metadata_df.withColumnRenamed("field.showcase", "phenocode")
phenotype_metadata_df.printSchema()

phenotypes_df = (
    filtered_genebass_df
    .sort("phenocode", "annotation")
    .groupby(
        # "annotation",
        "trait_type",
        "n_cases",
        "n_controls",
        "heritability",
        "phenocode",
        "pheno_sex",
        "coding",
        "modifier",
    )
    .agg(
        f.collect_set(f.col("annotation")).alias("annotation"),
        f.count(f.col("gene_id")).alias("num_significant_gene_associations"),
    )
    .persist()
    # .count()
    # .toPandas().astype({"n_controls": "Int32"})
)
phenotypes_df

filtered_phenotypes_df = (
    phenotypes_df
    .filter(
        (f.col("num_significant_gene_associations") > f.lit(1))
        & (f.col("annotation") != f.array(f.lit("synonymous")))
    )
    .withColumn(
        "total_num_samples", 
        f.col("n_cases") + f.when(f.col("n_controls").isNotNull(), f.col("n_controls")).otherwise(f.lit(0))
    )
)
filtered_phenotypes_df.toPandas()

joint_phenotypes_df = (
    filtered_phenotypes_df
    .join(phenotype_metadata_df.select("phenocode").distinct(), on="phenocode", how="inner")
    .filter(f.col("coding") == f.lit(""))
)
joint_phenotypes_df.toPandas().astype({"n_controls": "Int32"})

display(
    joint_phenotypes_df
    .join(phenotype_metadata_df, on="phenocode", how="inner")
    .toPandas().astype({"n_controls": "Int32"})
    .query("`col.name`.str.endswith('_0_0')")
    .loc[:, "col.name"].sort_values().values
)

phenocodes = joint_phenotypes_df.select("phenocode").distinct().toPandas()["phenocode"]
phenocodes.sort_values().values

len(phenocodes)

joint_phenotypes_df.toPandas().astype({"n_controls": "Int32"}).to_parquet(snakemake.output["filtered_phenotypes_pq"])


