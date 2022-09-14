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

# %%
from rep.notebook_init import init_ray
init_ray()

# %% {"tags": []}
import ray
from rep.notebook_init import init_spark_on_ray
spark = init_spark_on_ray(
    # executor_cores=128,
    # executor_memory_overhead=0.95,
    # configs={
    #     "spark.default.parallelism": int(ray.cluster_resources()["CPU"]),
    #     "spark.sql.shuffle.partitions": 2048,
    # }
)

# %%
spark

# %%
snakefile_path = os.getcwd() + "/../Snakefile"
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
        rule_name = 'associate_per_tissue',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0"
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))

# %%
n_genPCs = 30

# %%
age_col = 'age_when_attended_assessment_centre_f21003_0_0'
sex_col = 'sex_f31_0_0'
assert n_genPCs <= 40
genetic_principal_components_cols = [f"genetic_principal_components_f22009_0_{n}" for n in range(1, 1 + n_genPCs)]

# %%
covariate_cols = [
    age_col,
    sex_col,
    *genetic_principal_components_cols
#     'year_of_birth_f34_0_0',
]
covariate_cols

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %%
covariate_regression_formula = " + ".join([
    f'{phenotype_col} ~ 1',
    age_col,
    sex_col,
    f"{age_col} * {age_col}",
    f"{age_col} * {sex_col}",
    *genetic_principal_components_cols,
])
covariate_regression_formula

# %%
import re

m = re.search('.*_f(.+?)_.*', phenotype_col)
if m:
    phenocode = m.group(1)
else:
    raise ValueError("Cannot find phenocode!")
phenocode

# %%
# phenotype_column = (
#     phenotype_metadata_df
#     .query(f"phenocode == '{phenocode}'")
#     .loc[:, "col.name"]
# #     .iloc[:1]
#     .to_list()
# )
# phenotype_column

# %%
# age_fields = phenotype_metadata_df.loc[phenotype_metadata_df["field.tab"].str.contains('31')]
# age_fields

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
    *covariate_cols,
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
data_df = data_df.withColumnRenamed("eid", "sample_id")
# make sure that sample_id is a string
data_df = data_df.withColumn("sample_id", f.col("sample_id").cast(t.StringType()))

# broadcast it
data_df = f.broadcast(data_df.sort("sample_id"))

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
    covariate_regression_formula, 
    data = phenotype_df
).fit()

broadcast_phenotype_df = spark.sparkContext.broadcast(phenotype_df)
broadcast_restricted_model = spark.sparkContext.broadcast(restricted_model)


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
        phenotype_df = broadcast_phenotype_df.value
        
        # merge with phenotype df to make sure that we have all scores predicted
        data_df = phenotype_df.merge(pd_df, on=["sample_id"], how="left", ).fillna(0)
        
        with threadpool_limits(limits=1):
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
    dataframe=data_df.select("sample_id"), 
    groupby_columns=[],
    formula=covariate_regression_formula,
    restricted_formula=covariate_regression_formula,
    phenotype_df=broadcast_phenotype_df
).toPandas()

# %%
assert _test.iloc[0]["lr_df_diff"] == 0
display(_test)
del _test

# %% [markdown] {"tags": []}
# # read features

# %% [markdown]
# ## AbExp score

# %% [markdown]
# ### read AbExp predictions

# %%
abexp_sdf = (
    spark.read.parquet(snakemake.input["abexp_predictions"])
)
abexp_sdf.printSchema()

# %% {"tags": []}
abexp_sdf.select("subtissue").distinct().sort("subtissue").toPandas()

# %%
abexp_sdf.printSchema()

# %% {"tags": []}
test_df = (
    abexp_sdf
    .filter(f.col("gene") == f.lit('ENSG00000130173'))
    .toPandas()
)

# %%
test_df.query("gene_is_expressed")["subtissue"].unique()

# %% [markdown]
# ### perform regression

# %%
groupby_columns = ["gene", "subtissue"]

abexp_scores_sdf = (
    abexp_sdf
    # .filter(f.col("subtissue") == "Cells - Cultured fibroblasts")
    .filter(f.col("individual").isNotNull())
    .select(
        f.col("individual").alias("sample_id"),
        f.col("gene"),
        f.col("subtissue"),
        f.col("y_pred_proba").alias("AbExp_DNA"),
    )
)

abexp_scores_sdf.printSchema()

# %%
abexp_regression_results_sdf = regression(
    abexp_scores_sdf, 
    groupby_columns=["gene", "subtissue"], 
    formula=f"{covariate_regression_formula} + AbExp_DNA",
    restricted_formula=covariate_regression_formula,
    phenotype_df=broadcast_phenotype_df,
)
abexp_regression_results_sdf.printSchema()

# %% {"tags": []}
abexp_regression_results_sdf.write.parquet(snakemake.output["associations_pq"] + "/score_type=AbExp_DNA", mode="overwrite")

# %% [markdown] {"tags": []}
# ## pLoF counts

# %% [markdown] {"tags": []}
# ### read pLoF counts

# %%
plof_counts_sdf = (
    spark.read.parquet(snakemake.input["plof_counts"])
    .filter(f.col("sample_id").isNotNull())
    .select([
        f.col("sample_id").cast(t.StringType()), 
        f.col("gene_id").alias("gene"),
        "pLoF"
    ])
)
plof_counts_sdf.printSchema()

# %% [markdown] {"tags": []}
# ### perform regression

# %%
plof_regression_results_sdf = regression(
    plof_counts_sdf, 
    groupby_columns=["gene",], 
    formula=f"{covariate_regression_formula} + pLoF",
    restricted_formula=covariate_regression_formula,
    phenotype_df=broadcast_phenotype_df,
)
plof_regression_results_sdf.printSchema()

# %% {"tags": []}
plof_regression_results_sdf.write.parquet(snakemake.output["associations_pq"] + "/score_type=pLoF", mode="overwrite")

# %% [markdown]
# ## Joint regression

# %%
abexp_scores_sdf.printSchema()

# %%
plof_counts_sdf.printSchema()

# %%
joint_scores_df = (
    abexp_scores_sdf
    .join(plof_counts_sdf, on=["sample_id", "gene"], how="outer")
    .filter(f.col("sample_id").isNotNull())
    .filter(f.col("gene").isNotNull())
    .filter(f.col("subtissue").isNotNull())
    .fillna(0.)
)

# %%
joint_scores_df.printSchema()

# %%
joint_regression_results_sdf = regression(
    joint_scores_df, 
    groupby_columns=["gene", "subtissue",],
    formula=f"{covariate_regression_formula} + AbExp_DNA + pLoF",
    restricted_formula=covariate_regression_formula,
    phenotype_df=broadcast_phenotype_df,
)
joint_regression_results_sdf.printSchema()

# %% {"tags": []}
joint_regression_results_sdf.write.parquet(snakemake.output["associations_pq"] + "/score_type=joint", mode="overwrite")

# %% [markdown] {"tags": []}
# # compare pvalues

# %% [markdown]
# ## read protein-coding genes

# %%
protein_coding_genes_df = pd.read_parquet(snakemake.input["protein_coding_genes_pq"])#[["gene_id", "gene_name"]]
protein_coding_genes_df

# %% [markdown] {"tags": []}
# ## read association results

# %%
snakemake.output["associations_pq"] + "/score_type=AbExp_DNA"

# %%
abexp_regression_results_df = (
    spark.read.parquet(snakemake.output["associations_pq"] + "/score_type=AbExp_DNA")
    .sort("rsquared", ascending=False)
    .withColumn("score_pval", f.col("term_pvals")["AbExp_DNA"])
    .drop(
        "term_pvals",
        "params",
    )
    .toPandas()
)
abexp_regression_results_df
abexp_regression_results_df

# %% {"tags": []}
plof_regression_results_df = (
    spark.read.parquet(snakemake.output["associations_pq"] + "/score_type=pLoF")
    .sort("rsquared", ascending=False)
    .withColumn("score_pval", f.col("term_pvals")["pLoF"])
    .toPandas()
)

# %% {"tags": []}
# joint_regression_results_df = ( 
#     spark.read.parquet(snakemake.output["associations_pq"] + "/score_type=joint")
#     .filter(f.col("subtissue").isNotNull())
#     .sort("rsquared", ascending=False)
#     .toPandas()
# )

# %% {"tags": []}
# joint_regression_results_df_melted = ( 
#     spark.read.parquet(snakemake.output["associations_pq"] + "/score_type=joint")
#     .filter(f.col("subtissue").isNotNull())
#     .sort("rsquared", ascending=False)
#     .withColumn("term_pvals", f.expr("map_filter(`term_pvals`, (k, v) -> k IN ('pLoF', 'AbExp_DNA'))"))
#     .select("*", f.explode(f.col("term_pvals")).alias("score_type", "score_pval"))
#     .fillna(value=1, subset=["score_pval"])
#     .toPandas()
# )

# %% {"tags": []}
genebass_sdf = (
    spark.read.parquet("/s/project/rep/processed/genebass/results.parquet/")
    .filter(f.col("annotation") == f.lit("pLoF"))
    .filter(f.col("phenocode") == f.lit(phenocode))
)
genebass_sdf.printSchema()

# %% {"tags": []}
genebass_pd_df = genebass_sdf.toPandas()

# %%
genebass_pd_df

# %% [markdown]
# ### concat

# %%
adj_abexp_regression_results_df = (
    abexp_regression_results_df
    .loc[abexp_regression_results_df["gene"].isin(protein_coding_genes_df["gene_id"])]
#     .assign(score_type=("AbExp_DNA@" + adj_abexp_regression_results_df["subtissue"]).astype("string[pyarrow]"))
#     .drop(columns=["subtissue"])
    .assign(score_type="AbExp-DNA")
    .astype({"score_type": "string[pyarrow]"})
)

print(f"Correcting for {adj_abexp_regression_results_df.shape[0]} association tests...")
adj_abexp_regression_results_df = adj_abexp_regression_results_df.assign(padj=np.fmin(adj_abexp_regression_results_df["lr_pval"] * adj_abexp_regression_results_df.shape[0], 1))

adj_abexp_regression_results_df

# %%
(
    pn.ggplot(adj_abexp_regression_results_df.query("padj < 0.05"), pn.aes(x="subtissue"))
    + pn.geom_bar()
    + pn.theme(axis_text_x=pn.element_text(rotation=90))
    + pn.ggtitle("Number of significantly associating genes per subtissue")
)

# %%
# keep only subtissue with most significant padj for each gene
adj_abexp_regression_results_df = adj_abexp_regression_results_df.groupby("gene").first().sort_values("padj").reset_index()

# %%
adj_abexp_regression_results_df["padj"].quantile(q=[0.5, 0.05, 0.001])

# %%
(adj_abexp_regression_results_df["padj"] < 0.05).sum()

# %%
adj_plof_regression_results_df = (
    plof_regression_results_df
    .assign(score_type="pLoF")
    .astype({"score_type": "string[pyarrow]"})
)
adj_plof_regression_results_df["padj"] = np.fmin(adj_plof_regression_results_df["lr_pval"] * adj_plof_regression_results_df.shape[0], 1)
adj_plof_regression_results_df

# %%
adj_plof_regression_results_df["padj"].quantile(q=[0.5, 0.05])

# %%
(adj_plof_regression_results_df["padj"] < 0.05).sum()

# %%
adj_genebass_regression_results_df = (
    genebass_pd_df
    .assign(score_type="pLoF")
    .rename(columns={
        "gene_id": "gene",
        "Pvalue": "lr_pval",
    })
    .assign(score_type="Genebass")
    .astype({"score_type": "string[pyarrow]"})
)
adj_genebass_regression_results_df["padj"] = np.fmin(adj_genebass_regression_results_df["lr_pval"] * adj_genebass_regression_results_df.shape[0], 1)
adj_genebass_regression_results_df

# %%
adj_genebass_regression_results_df["padj"].quantile(q=[0.5, 0.05])

# %%
(adj_genebass_regression_results_df["padj"] < 0.05).sum()

# %%
# adj_joint_regression_results_df = (
#     joint_regression_results_df
#     .loc[joint_regression_results_df["gene"].isin(protein_coding_genes_df["gene_id"])]
#     .assign(score_type=("AbExp-DNA@" + joint_regression_results_df["subtissue"] + " + pLoF").astype("string[pyarrow]"))
#     .drop(columns=["subtissue"])
# )
# adj_joint_regression_results_df["padj"] = np.fmin(adj_joint_regression_results_df["lr_pval"] * adj_joint_regression_results_df.shape[0], 1)
# adj_joint_regression_results_df

# %%
# adj_joint_regression_results_df["padj"].quantile(q=[0.5, 0.05, 0.001])

# %%
# (adj_joint_regression_results_df["padj"] < 0.05).sum()

# %%
regression_results_df = pd.concat([
    adj_abexp_regression_results_df,
    adj_plof_regression_results_df,
    adj_genebass_regression_results_df,
    # adj_joint_regression_results_df,
], join="inner")
regression_results_df = regression_results_df.set_index(["gene", "score_type"])["padj"].unstack("score_type")
# regression_results_df = regression_results_df.fillna(1)
regression_results_df = regression_results_df.sort_values("Genebass")

# filiter for Genebass genes
regression_results_df = regression_results_df.dropna(subset=["Genebass"])
regression_results_df

# %% [markdown]
# ## Plot

# %%
x="Genebass"
y="AbExp-DNA"

cutoff = 0.05

counts_x = (regression_results_df[x] < cutoff).sum()
counts_y = (regression_results_df[y] < cutoff).sum()

min_pval = min(
    regression_results_df[x].min(),
    regression_results_df[y].min()
)

(
    pn.ggplot(regression_results_df.fillna(1), pn.aes(x=x, y=y))
    + pn.geom_abline(slope=1, color="black", linetype="dashed")
    + pn.geom_hline(yintercept=0.05, color="red", linetype="dashed")
    + pn.geom_vline(xintercept=0.05, color="red", linetype="dashed")
    # + pn.geom_rect(
    #     xmax=1,
    #     xmin=0.05,
    #     ymax=1,
    #     ymin=0.05,
    #     color="red", 
    #     linetype="dashed"
    # )
    + pn.geom_point()
    #+ pn.coord_fixed()
    + pn.scale_x_log10(limits=(min_pval, 1))
    + pn.scale_y_log10(limits=(min_pval, 1))
    + pn.labs(
        title=f"Association between genes and '{phenotype_col}'\n(p-values, alpha={cutoff})",
        x=f"{x}\n(n_signif={counts_x})",
        y=f"{y}\n(n_signif={counts_y})",
    )
    # + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
)

# %%
with pd.option_context('display.float_format', '{:,.2g}'.format):
    display(
        (
            regression_results_df
            .loc[(regression_results_df[x] < cutoff) | (regression_results_df[y] < cutoff), [x, y]]
            .merge(protein_coding_genes_df.rename(columns={"gene_id": "gene"}), on="gene", how="left")
            .set_index(["gene", "gene_name"])
            .sort_values(x)
        )
    )

# %%
x="pLoF"
y="AbExp-DNA"

cutoff = 0.05

counts_x = (regression_results_df[x] < cutoff).sum()
counts_y = (regression_results_df[y] < cutoff).sum()

min_pval = min(
    regression_results_df[x].min(),
    regression_results_df[y].min()
)

(
    pn.ggplot(regression_results_df.fillna(1), pn.aes(x=x, y=y))
    + pn.geom_abline(slope=1, color="black", linetype="dashed")
    + pn.geom_hline(yintercept=0.05, color="red", linetype="dashed")
    + pn.geom_vline(xintercept=0.05, color="red", linetype="dashed")
    # + pn.geom_rect(
    #     xmax=1,
    #     xmin=0.05,
    #     ymax=1,
    #     ymin=0.05,
    #     color="red", 
    #     linetype="dashed"
    # )
    + pn.geom_point()
    #+ pn.coord_fixed()
    + pn.scale_x_log10(limits=(min_pval, 1))
    + pn.scale_y_log10(limits=(min_pval, 1))
    + pn.labs(
        title=f"Association between genes and '{phenotype_col}'\n(p-values, alpha={cutoff})",
        x=f"{x} p-value\n(n_signif={counts_x})",
        y=f"{y} p-value\n(n_signif={counts_y})",
    )
    # + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
)

# %%
abexp_regression_results_df.subtissue.unique().shape

# %%
with pd.option_context('display.float_format', '{:,.2g}'.format):
    display(
        (
            regression_results_df
            .loc[(regression_results_df[x] < cutoff) | (regression_results_df[y] < cutoff), [x, y]]
            .merge(protein_coding_genes_df.rename(columns={"gene_id": "gene"}), on="gene", how="left")
            .set_index(["gene", "gene_name"])
            .sort_values(x)
            # .fillna(1)
        )
    )

# %%

# %%

# %%

# %%

