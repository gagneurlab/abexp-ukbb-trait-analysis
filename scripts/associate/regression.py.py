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
from IPython.display import display

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

# %%
import plotnine as pn
# import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

# %%
# import os
# # os.environ["RAY_ADDRESS"] = os.environ.get("RAY_ADDRESS", 'ray://192.168.16.30:10001')
# os.environ["RAY_ADDRESS"] = 'ray://192.168.16.28:10001'
# os.environ["RAY_ADDRESS"]

# %% {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark()

# %% [raw]
# from rep.notebook_init import init_ray
# init_ray()

# %% [raw] {"tags": []}
# import ray
# from rep.notebook_init import init_spark_on_ray
# spark = init_spark_on_ray(
#     # executor_cores=128,
#     executor_memory_overhead=0.9,
#     # configs={
#     #     "spark.default.parallelism": int(ray.cluster_resources()["CPU"]),
#     #     "spark.sql.shuffle.partitions": 2048,
#     # }
# )

# %%
spark

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
        rule_name = 'associate__regression',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            "feature_set": "LOFTEE_pLoF",
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))

# %% [markdown]
# # Load configuration

# %%
with open(snakemake.input["featureset_config"], "r") as fd:
    config = yaml.safe_load(fd)

# %%
print(json.dumps(config, indent=2, default=str))

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %%
phenocode = config["phenocode"]
phenocode

# %%
restricted_formula = config["restricted_formula"]
print(restricted_formula)

# %% [markdown]
# # Read metadata

# %%
phenotype_metadata_df = (
    spark.read.parquet(snakemake.input["phenotype_metadata_pq"])
    .withColumnRenamed("field.showcase", "phenocode")
    .toPandas()
)
phenotype_metadata_df

# %%
phenotype_metadata_subset = phenotype_metadata_df.set_index("col.name").loc[[
    *config["covariate_cols"],
#     *phenotype_column,
    phenotype_col,
]]
phenotype_metadata_subset

# %%
phenotype_metadata_subset.reset_index().to_parquet(snakemake.output["metadata_pq"], index=False)

# %%
phenotype_metadata_subset.reset_index().to_csv(snakemake.output["metadata_tsv"], index=False, header=True, sep="\t")

# %% [markdown]
# # Load phenotypes

# %%
data_dfs = []
for data_path, group_df in phenotype_metadata_subset.groupby("data_path"):
    data_df = (
        spark.read.parquet(data_path)
        .select(
            "eid",
            *group_df.index.to_list()
        )
    )
    data_dfs.append(data_df)

# %%
len(data_dfs)

# %%
import functools
data_df = functools.reduce(lambda df1, df2: df1.join(df2, on="eid", how="outer"), data_dfs)

# remove NA values from the phenotype column
data_df = data_df.filter(f.col(phenotype_col).isNotNull())
# rename eid col
data_df = data_df.withColumnRenamed("eid", "individual")
# make sure that sample_id is a string
data_df = data_df.withColumn("individual", f.col("individual").cast(t.StringType()))

# broadcast it
data_df = f.broadcast(data_df.sort("individual"))

data_df.printSchema()

# %% {"tags": []}
spark._jvm.System.gc()

# %% {"tags": []}
phenotype_df = data_df.toPandas()

# %%
phenotype_df

# %% {"tags": []}
import statsmodels.formula.api as smf
from threadpoolctl import threadpool_limits
from typing import List, Union

restricted_model = smf.ols(
    restricted_formula, 
    data = phenotype_df
).fit()

broadcast_phenotype_df = spark.sparkContext.broadcast(phenotype_df)
# broadcast_restricted_model = spark.sparkContext.broadcast(restricted_model)

# %% {"tags": []}
def regression(
    dataframe: pyspark.sql.DataFrame,
    groupby_columns: List[str], 
    formula: str,
    restricted_formula: str,
    phenotype_df: Union[pyspark.Broadcast, pd.DataFrame]
):
    # make sure to use broadcasted variable
    if isinstance(phenotype_df, pyspark.Broadcast):
        broadcast_phenotype_df = phenotype_df
    elif isinstance(phenotype_df, pd.DataFrame):
        broadcast_phenotype_df = spark.sparkContext.broadcast(phenotype_df)
    else:
        raise ValueError("'phenotype_df' has to be of type pd.DataFrame or pyspark.Broadcast")
    
#     if isinstance(restricted_model, pyspark.Broadcast):
#         broadcast_restricted_model = restricted_model
#     else:
#         broadcast_restricted_model = spark.sparkContext.broadcast(restricted_model)
    
    def fit(pd_df):
        with threadpool_limits(limits=1):
            phenotype_df = broadcast_phenotype_df.value

            # merge with phenotype df to make sure that we have all scores predicted
            data_df = phenotype_df.merge(pd_df, on=["individual"], how="left", ).fillna(0)
        
            # restricted_model = broadcast_restricted_model.value
            restricted_model = smf.ols(
                restricted_formula,
                data = data_df
            ).fit()
            
#             randomized_model = smf.ols(
#                 restricted_formula,
#                 data = data_df
#             ).fit()
            
            model = smf.ols(
                formula,
                data = data_df
            ).fit()
            
            lr_stat, lr_pval, lr_df_diff = model.compare_lr_test(restricted_model)

            return (
                pd_df
                .iloc[:1].loc[:, groupby_columns]
                .copy()
                .assign(**{
                    "n_samples": [data_df.shape[0]],
                    "term_pvals": [model.pvalues.to_dict()], 
                    "params": [model.params.to_dict()], 
                    "loglikelihood": [model.llf],
                    "rsquared_restricted": [restricted_model.rsquared],
                    "rsquared": [model.rsquared],
                    "lr_stat": [lr_stat],
                    "lr_pval": [lr_pval],
                    "lr_df_diff": [lr_df_diff],
                })
            )
    
    return dataframe.groupby(groupby_columns).applyInPandas(
        func=fit,
        schema=t.StructType([
            *[dataframe.schema[k] for k in groupby_columns],
            t.StructField("n_samples", t.LongType()),
            t.StructField("term_pvals", t.MapType(t.StringType(), t.DoubleType())),
            t.StructField("params", t.MapType(t.StringType(), t.DoubleType())),
            t.StructField("loglikelihood", t.DoubleType()),
            t.StructField("rsquared_restricted", t.DoubleType()),
            t.StructField("rsquared", t.DoubleType()),
            t.StructField("lr_stat", t.DoubleType()),
            t.StructField("lr_pval", t.DoubleType()),
            t.StructField("lr_df_diff", t.DoubleType()),
        ])
    )

_test = regression(
    dataframe=data_df.select("individual"), 
    groupby_columns=[],
    formula=restricted_formula,
    restricted_formula=restricted_formula,
    phenotype_df=broadcast_phenotype_df
).toPandas()

# %%
assert _test.iloc[0]["lr_df_diff"] == 0
display(_test)
del _test

# %% [markdown] {"tags": []}
# # read features

# %%
groupby_columns = config["groupby_columns"]
groupby_columns

# %%
config["feature_sets"]

# %%
config["variables"]

# %% [markdown]
# ## protein-coding genes

# %%
protein_coding_genes_df = (
    spark.read.parquet(snakemake.input["protein_coding_genes_pq"])
    .withColumnRenamed("gene_id", "gene")
)
protein_coding_genes_df.printSchema()

# %% [markdown]
# ## feature dfs

# %%
feature_dfs = {}
for feature_name, path in config["feature_sets"].items():
    feature_dfs[feature_name] = spark.read.parquet(path + "/data.parquet")

# %%
len(feature_dfs)

# %%
from rep.data import join_featuresets

# %%
fill_values = config.get('fill_values')
if fill_values is not None and len(fill_values) > 0:
    # quote the column names with backticks
    fill_values = {"`" + k + "`": v for k, v in fill_values.items()}

features_df = join_featuresets(
    dataframes=feature_dfs,
    variables=config["variables"],
    index_cols=["individual", *groupby_columns],
    fill_values=fill_values,
    join="left",
    initial_df=protein_coding_genes_df
).filter(f.col("individual").isNotNull())
features_df.printSchema()

# %%
renamed_features_df = features_df
features = []

for c in features_df.columns:
    if c.startswith("feature."):
        # rename columns because spark is stupid and fails with selecting in a pandas udf
        new_name = c[8:].replace(".", "_").replace("@", "__")
        renamed_features_df = renamed_features_df.withColumnRenamed(c, new_name)
        
        features.append(new_name)

renamed_features_df.printSchema()

# %%
features

# %% [markdown]
# # perform regression

# %%
# scores_sdf = (
#     features_df
#     # .filter(f.col("subtissue") == "Cells - Cultured fibroblasts")
#     .filter(f.col("individual").isNotNull())
#     .select(
#         "individual",
#         *groupby_columns,
#         *[f"`{c}`" for c in features],
#     )
# )

# scores_sdf.printSchema()

# %%
full_formula = "\n + ".join([
    restricted_formula,
    *features,
])
print(full_formula)

# %%
regression_results_sdf = regression(
    renamed_features_df, 
    groupby_columns=groupby_columns, 
    formula=full_formula,
    restricted_formula=restricted_formula,
    phenotype_df=broadcast_phenotype_df,
)
regression_results_sdf.printSchema()

# %% {"tags": []}
regression_results_sdf.write.parquet(snakemake.output["associations_pq"], mode="overwrite")

# %% [raw] {"tags": []}
# # read association results

# %% [raw]
# snakemake.output["associations_pq"]

# %% [raw]
# regression_results_df = (
#     spark.read.parquet(snakemake.output["associations_pq"])
#     .sort("rsquared", ascending=False)
#     # .withColumn("score_pval", f.col("term_pvals")["AbExp_DNA"])
#     .drop(
#         "term_pvals",
#         "params",
#     )
#     .toPandas()
# )
# regression_results_df

# %%

