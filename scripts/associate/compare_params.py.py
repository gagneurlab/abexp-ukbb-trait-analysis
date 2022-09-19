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

# %% {"tags": []}
from rep.notebook_init import init_spark
spark = init_spark()

# %%
spark

# %%
# import os
# # os.environ["RAY_ADDRESS"] = os.environ.get("RAY_ADDRESS", 'ray://192.168.16.30:10001')
# os.environ["RAY_ADDRESS"] = 'ray://192.168.16.28:10001'
# os.environ["RAY_ADDRESS"]

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %%
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__compare_params',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            # "feature_set": "LOFTEE_pLoF",
            "feature_set": "AbExp_all_tissues",
            "covariates": "sex_age_genPC",
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))

# %%
if "plot_dpi" in snakemake.params:
    DPI = snakemake.params["plot_dpi"]
else:
    DPI=450

# %% [markdown] {"tags": []}
# # read features

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

# %% [markdown]
# ## read protein-coding genes

# %%
protein_coding_genes_df = spark.read.parquet(snakemake.input["protein_coding_genes_pq"])
protein_coding_genes_df.limit(10).toPandas()

# %% [markdown] {"tags": []}
# ## read association results

# %%
snakemake.input["associations_pq"]

# %%
regression_results_df = spark.read.parquet(snakemake.input["associations_pq"])

# %%
correction_value = regression_results_df.count()
print(f"Corrected for {correction_value} association tests...")

# %%
regression_results_df = (
    regression_results_df
    .sort("rsquared", reverse=True)
    .withColumn("padj", f.array_min(f.array(
        f.col("lr_pval") * f.lit(correction_value),
        f.lit(1.0),
    )))
    .withColumn("rsquared_diff", f.col("rsquared") - f.col("rsquared_restricted"))
    # .collect()
    # .to_pandas()
)

# %%
regression_results_df.printSchema()

# %%

f.udf(

# %%
from urllib.parse import unquote

spark_unquote = f.udf(unquote, t.StringType())

# %%
regression_params_df = (
    regression_results_df
    .select(
        "*",
        f.explode(f.col("term_pvals")).alias("term", "term_pval")
    )
    .withColumn("term", spark_unquote(f.regexp_extract(f.col("term"), r'''(Q\(['"])?([^"']*)''', 2)))
    .filter(f.col("term").contains("AbExp"))
    .withColumn("term", f.regexp_replace(f.col("term"), r'''AbExp@AbExp_''', ""))
    .drop(
        "term_pvals",
        "params",
    )
)
regression_params_df.printSchema()

# %%
terms = regression_params_df.select("term").distinct().toPandas()["term"].sort_values()
terms

# %%
regression_params_pd_df = regression_params_df.toPandas()
regression_params_pd_df

# %%

# %%

# %% [markdown]
# ## Plot

# %% [markdown]
# ## full plot

# %%
plot_df = regression_params_pd_df
plot_df = (
    plot_df
    .assign(minus_log10_pval=-np.log10(plot_df["term_pval"].values))
    .query("padj < 0.05")
)
plot_df

# %%
box_order = plot_df.groupby("term")["minus_log10_pval"].median().sort_values()
box_order

# %%
plot = (
    pn.ggplot(plot_df, pn.aes(
        x="term",
        y="minus_log10_pval",
    ))
    + pn.geom_boxplot()
    + pn.theme(axis_text_x=pn.element_text(rotation=90))
    + pn.scale_x_discrete(limits=box_order.index)
    + pn.labs(
        title=f"Distribution of term p-values for significantly associating genes\n(trait: '{phenotype_col}')",
        #x="",
        y="-log10(p-value)",
    )
)

# %%
display(plot)

# %%
path = snakemake.params["output_basedir"] + "/term_pvalue_boxplot"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %%

