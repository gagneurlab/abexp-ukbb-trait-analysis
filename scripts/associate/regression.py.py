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
            # "feature_set": "LOFTEE_pLoF",
            # "feature_set": "AbExp_pivot",
            "feature_set": "LOFTEE_pLoF",
            # "covariates": "sex+age+genPC+CLMP",
            "covariates": "sex+age+genPC",
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
phenocode = config["covariates"]["phenocode"]
phenocode

# %%
restricted_formula = config["covariates"]["restricted_formula"]
print(restricted_formula)

# %% [markdown]
# # Read covariates

# %%
covariates_df = pd.read_parquet(snakemake.input["covariates_pq"]).fillna(0)
covariates_df

# %%
display(
    "Size of 'covariates_df': %.3fGb" % (covariates_df.memory_usage(deep=True).sum() / 1024**3)
)

# %% [markdown] {"tags": []}
# ## clumping

# %%
clumping_variants_df = pd.read_parquet(snakemake.input["clumping_variants_pq"])
clumping_variants_df


# %%
def get_variants_by_gene(clumping_variants_df, gene_id):
    included_vars = clumping_variants_df.query(f"gene == '{gene_id}'")["variant"].values
    return included_vars


# %%
# get_variants_by_gene(clumping_variants_df, gene_id="ENSG00000084674")

# %%
import re

def format_formula(formula, keys, add_clumping=True, clumping_variants_df=clumping_variants_df):
    if add_clumping:
        if not "gene" in keys:
            raise ValueError(f"missing gene in keys: '{keys}'!")

        gene_id = keys["gene"]

        variants = get_variants_by_gene(clumping_variants_df=clumping_variants_df, gene_id=gene_id)

        if len(variants) > 0:
            formula = "\n + ".join([
                formula,
                *[f"Q('{c}')" for c in variants]
            ])
    
    return formula


# %%
test_formula = format_formula(
    formula=restricted_formula,
    clumping_variants_df=clumping_variants_df,
    add_clumping=True,
    keys={
        "gene": "ENSG00000084674",
    }
)

# %%
# print(test_formula)

# %%
test_formula_2 = format_formula(
    formula=restricted_formula,
    clumping_variants_df=clumping_variants_df,
    add_clumping=True,
    keys={
        "gene": "",
    }
)

# %%
# print(test_formula_2)

# %%
import patsy
import re

def get_variables_from_formula(formula, lhs=True, rhs=True):
    model_desc = patsy.ModelDesc.from_formula(formula)

    covariate_cols = [
        *([factor.name() for term in model_desc.lhs_termlist for factor in term.factors] if lhs else []),
        *([factor.name() for term in model_desc.rhs_termlist for factor in term.factors] if rhs else []),
    ]
    # deduplicate
    covariate_cols = list(dict.fromkeys(covariate_cols))
    
    # replace Q($i) -> $i
    regex = re.compile(r"""Q\(['"](.*)['"]\)""")
    covariate_cols = [regex.sub(r"\1", c) for c in covariate_cols]
    
    return covariate_cols


# %%
# get_variables_from_formula(test_formula)

# %% {"tags": []}
broadcast_covariates_df = spark.sparkContext.broadcast(covariates_df)
broadcast_clumping_variants_df = spark.sparkContext.broadcast(clumping_variants_df)

# %% {"tags": []}
# restricted_model = smf.ols(
#     restricted_formula,
#     data = (
#         covariates_df
#         .assign(sex_f31_0_0=covariates_df["sex_f31_0_0"].where(np.random.randint(0, 2, size=covariates_df.shape[0], dtype=bool))),
#     )
# ).fit()
# broadcast_restricted_model = spark.sparkContext.broadcast(restricted_model)

# %% {"tags": []}
import statsmodels.formula.api as smf

from threadpoolctl import threadpool_limits
from typing import List, Union

def regression(
    dataframe: pyspark.sql.DataFrame,
    groupby_columns: Union[str, List[str]], 
    full_formula: str,
    restricted_formula: str,
    covariates_df: Union[pyspark.Broadcast, pd.DataFrame],
    clumping_variants_df: Union[pyspark.Broadcast, pd.DataFrame] = None,
    add_clumping: bool = config["covariates"]["add_clumping"],
):
    if isinstance(groupby_columns, str):
        groupby_columns = [groupby_columns]
    
    # make sure to use broadcasted variable
    if isinstance(covariates_df, pyspark.Broadcast):
        broadcast_covariates_df = covariates_df
    elif isinstance(covariates_df, pd.DataFrame):
        broadcast_covariates_df = spark.sparkContext.broadcast(covariates_df)
    else:
        raise ValueError("'covariates_df' has to be of type pd.DataFrame or pyspark.Broadcast")
    
    if add_clumping:
        if isinstance(clumping_variants_df, pyspark.Broadcast):
            broadcast_clumping_variants_df = clumping_variants_df
        elif isinstance(clumping_variants_df, pd.DataFrame):
            broadcast_clumping_variants_df = spark.sparkContext.broadcast(clumping_variants_df)
        else:
            raise ValueError("'clumping_variants_df' has to be of type pd.DataFrame or pyspark.Broadcast")
    else:
        broadcast_clumping_variants_df = None
    
#     if isinstance(restricted_model, pyspark.Broadcast):
#         broadcast_restricted_model = restricted_model
#     else:
#         broadcast_restricted_model = spark.sparkContext.broadcast(restricted_model)
    
    def fit(pd_df):
        with threadpool_limits(limits=1):
            keys=pd_df.loc[:, groupby_columns].iloc[0].to_dict()
            
            # ## unpack all broadcast variables
            covariates_df = broadcast_covariates_df.value
            
            # unpack only if needed
            if add_clumping:
                clumping_variants_df = broadcast_clumping_variants_df.value
            else:
                clumping_variants_df = None
            # ## unpacking done

            formatted_restricted_formula = format_formula(
                formula=restricted_formula,
                keys=keys,
                add_clumping=add_clumping,
                clumping_variants_df=clumping_variants_df,
            )
            formatted_full_formula = format_formula(
                formula=full_formula,
                keys=keys,
                add_clumping=add_clumping,
                clumping_variants_df=clumping_variants_df,
            )
            
            # get necessary columns
            restricted_variables = get_variables_from_formula(formatted_restricted_formula)
            full_variables = get_variables_from_formula(formatted_full_formula)
            necessary_columns = {
                *restricted_variables,
                *full_variables,
                "individual",
            }
            # filter covariates for necessary columns
            covariates_df = covariates_df.loc[:, [c for c in covariates_df.columns if c in necessary_columns]]
            
            # merge with phenotype df to make sure that we have all scores predicted
            data_df = covariates_df.merge(pd_df, on=["individual"], how="left")
            # fill missing values
            data_df = data_df.fillna({
                c: 0 for c in pd_df.columns
            })
        
            # restricted_model = broadcast_restricted_model.value
            restricted_model = smf.ols(
                formatted_restricted_formula,
                data = data_df
            ).fit()
            
#             randomized_model = smf.ols(
#                 restricted_formula,
#                 data = data_df
#             ).fit()
            
            model = smf.ols(
                formatted_full_formula,
                data = data_df
            ).fit()
            
            lr_stat, lr_pval, lr_df_diff = model.compare_lr_test(restricted_model)

            return (
                pd_df
                .loc[:, groupby_columns]
                .iloc[:1]
                .copy()
                .assign(**{
                    "n_observations": [int(model.nobs)],
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
            t.StructField("n_observations", t.LongType()),
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
    dataframe=spark.createDataFrame(covariates_df[["individual"]].assign(gene="ENSG00000084674")), 
    groupby_columns=["gene"],
    full_formula=restricted_formula,
    restricted_formula=restricted_formula,
    covariates_df=broadcast_covariates_df,
    clumping_variants_df=broadcast_clumping_variants_df,
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
        new_name = (
            c[8:]
            .replace(".", "_")
            # .replace("@", "__")
        )
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
    *[f"Q('{c}')" for c in features],
])
print(full_formula)

# %%
regression_results_sdf = regression(
    renamed_features_df, 
    groupby_columns=groupby_columns, 
    full_formula=full_formula,
    restricted_formula=restricted_formula,
    covariates_df=broadcast_covariates_df,
    clumping_variants_df=broadcast_clumping_variants_df,
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

