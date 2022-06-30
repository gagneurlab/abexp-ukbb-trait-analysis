# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.2
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

# %% {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark()

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
        rule_name = 'associate',
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
data_df.printSchema()

# %%
# remove NA values from the phenotype column
data_df = data_df.filter(f.col(phenotype_col).isNotNull())
# rename eid col
data_df = data_df.withColumnRenamed("eid", "sample_id")
# make sure that sample_id is a string
data_df = data_df.withColumn("sample_id", f.col("sample_id").cast(t.StringType()))

# %%
data_pd_df = data_df.toPandas().astype({"sample_id": "string[pyarrow]"})
data_pd_df

# %%
data_pd_df.to_parquet(snakemake.output["data_pq"], index=False)

# %% [markdown]
# # regress out standard covariates

# %%
import statsmodels.formula.api as smf

model = smf.ols(
    covariate_regression_formula, 
    data = data_pd_df
).fit()

# %%
print(model.params)

# %%
model.rsquared

# %%
data_pd_df["residuals"] = data_pd_df[phenotype_col] - model.predict(data_pd_df)
data_pd_df

# %%
data_pd_df["residuals"].hist()

# %% [markdown]
# # read features

# %% [markdown]
# ## read AbExp predictions

# %%
abexp_sdf = (
    spark.read.parquet('/s/project/rep/processed/ukbb_wes_200k/genebass_abexp_pred.parquet')
)
abexp_sdf.printSchema()

# %%
abexp_sdf.select("subtissue").distinct().sort("subtissue").toPandas()

# %%
abexp_df = (
    abexp_sdf
    .filter(f.col("subtissue") == "Cells - Cultured fibroblasts")
    .withColumnRenamed("individual", "sample_id")
    .withColumnRenamed("gene", "gene_id")
    .sort("gene_id", "sample_id")
    .toPandas()
)
abexp_df

# %%
len(abexp_df.sample_id.unique())

# %%
abexp_df.y_pred_proba.max()

# %%
abexp_df.y_pred_proba.min()

# %%
abexp_df.y_pred_proba.hist(bins=50, log=True)

# %%
abexp_scores_df = (
    abexp_df
    .rename(columns={"y_pred_proba": "score"})
    .assign(score_type="AbExp-DNA")
    .loc[:, ["sample_id", "gene_id", "score", "score_type"]]
)
abexp_scores_df

# %% [markdown]
# ## read pLoF counts

# %%
plof_counts_df = (
    pd.read_parquet(snakemake.input["plof_counts"])
    .rename(columns={"pLoF": "score"})
    .assign(score_type="pLoF")
    .loc[:, ["sample_id", "gene_id", "score", "score_type"]]
)
plof_counts_df

# %%
plof_counts_df.query("score > 0")

# %%
plof_counts_df.score.hist(bins=50, log=True)

# %% [markdown]
# ## concat all scores

# %%
scores_df = pd.concat([
    abexp_scores_df,
    plof_counts_df,
], axis=0)
scores_df = scores_df.astype({
    "sample_id": "string[pyarrow]",
    "score_type": "string[pyarrow]",
})
scores_df["rank"] = scores_df.groupby("score_type")["score"].rank(method="max", ascending=False).astype("int32")
scores_df

# %%
target_rank = scores_df.query(f"score_type == 'pLoF'").query("score > 0")["rank"].max()
target_rank

# %%
scores_df["is_significant"] = scores_df["rank"] <= target_rank
scores_df

# %%
print("Number of significants:")
print(scores_df.groupby("score_type")["is_significant"].sum().to_dict())

# %% [markdown]
# # read genebass

# %%
genebass_df = (
    spark.read.parquet(snakemake.input["genebass_pq"])
    .filter(f.col("annotation") == f.lit("pLoF"))
    .filter(f.col("phenocode") == f.lit(phenocode))
)
genebass_df.printSchema()

# %%
genebass_pd_df = genebass_df.toPandas()
genebass_pd_df

# %%
# phenotype_metadata_df.merge(genebass_pd_df, on="phenocode", how="inner")

# %%
plot_df = (
    data_pd_df[["sample_id", "residuals"]]
    .merge(
        genebass_pd_df[["annotation", "phenocode", "gene_id", "gene_symbol", "total_variants"]].merge(scores_df, on="gene_id", how="inner"),
        on="sample_id"
    )
)
plot_df

# %% pn.facet_wrap(["gene_id", "score_type"])
# plot = (
#     pn.ggplot(plot_df, pn.aes(x="score", y="residuals"))
#     + pn.stat_bin2d()
#     + pn.facet_wrap(["gene_id", "score_type"])
# )
# plot


# %%
plot = (
    pn.ggplot(plot_df, pn.aes(x="score", y="residuals"))
    + pn.stat_bin2d()
    + pn.facet_grid("gene_symbol ~ score_type", scales="free_x")
    # + pn.coord_fixed(ratio = 1/2)
    + pn.theme(aspect_ratio=1/2)
)
plot

# %%
plot = (
    pn.ggplot(plot_df, pn.aes(x="score", y="residuals"))
    + pn.stat_bin2d()
    + pn.facet_grid("gene_symbol ~ score_type", scales="free_x")
    + pn.coord_fixed(ratio = 1/4)
)
plot

# %%

# %% {"incorrectly_encoded_metadata": "pn.theme(aspect_ratio=1/2)"}
plot = (
    pn.ggplot(plot_df, pn.aes(x="score_type", y="residuals", color="is_significant"))
    + pn.geom_boxplot()
    + pn.theme(axis_text_x=pn.element_text(rotation=90, hjust=1))
    + pn.labs(
    #         x="gene id",
    #         y="residuals",
        title=f"{phenotype_col}"
    )
    + pn.facet_grid("~ gene_symbol")
    + pn.theme(figure_size=(8, 6))
)
plot

# %%
plot 

# %% {"incorrectly_encoded_metadata": "pn.theme(aspect_ratio=1/2)"}
plot = (
    pn.ggplot(plot_df, pn.aes(x="score_type", y="residuals", color="is_significant"))
    + pn.geom_boxplot()
#     + pn.theme(axis_text_x=pn.element_text(rotation=90, hjust=1))
    + pn.labs(
    #         x="gene id",
    #         y="residuals",
        title=f"{phenotype_col}"
    )
    + pn.facet_grid("gene_symbol ~")
    + pn.theme(figure_size=(6, 8)) 
    + pn.coord_flip()
)
plot

# %%

