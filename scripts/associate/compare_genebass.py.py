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
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__compare_genebass',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            "feature_set": "LOFTEE_pLoF",
            "covariates": "sex+age+genPC",
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
protein_coding_genes_df = pd.read_parquet(snakemake.input["protein_coding_genes_pq"])#[["gene_id", "gene_name"]]
protein_coding_genes_df

# %% [markdown] {"tags": []}
# ## read association results

# %%
snakemake.input["associations_pq"]

# %%
regression_results_df = (
    spark.read.parquet(snakemake.input["associations_pq"])
    .sort("rsquared", ascending=False)
    .drop(
        "term_pvals",
        "params",
    )
    .toPandas()
)

print(f"Correcting for {regression_results_df.shape[0]} association tests...")
regression_results_df = regression_results_df.assign(padj=np.fmin(regression_results_df["lr_pval"] * regression_results_df.shape[0], 1))

display(regression_results_df)

# %%
regression_results_df["padj"].quantile(q=[0.5, 0.05, 0.001])

# %%
(regression_results_df["padj"] < 0.05).sum()

# %%
# regression_results_df[regression_results_df["padj"] < 0.05]

# %% {"tags": []}
genebass_sdf = (
    spark.read.parquet(snakemake.input["genebass_pq"])
    .filter(f.col("annotation") == f.lit("pLoF"))
    .filter(f.col("phenocode") == f.lit(phenocode))
)
genebass_sdf.printSchema()

# %% {"tags": []}
genebass_pd_df = genebass_sdf.toPandas()

# %%
genebass_pd_df

# %%
adj_genebass_regression_results_df = (
    genebass_pd_df
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
combined_regression_results_df = pd.concat([
    (
        regression_results_df
        .assign(score_type=snakemake.wildcards["feature_set"])
        .astype({"score_type": "string[pyarrow]"})
    ),
    adj_genebass_regression_results_df,
], join="inner")

combined_regression_results_df = (
    combined_regression_results_df
    .set_index(["gene", "score_type"])["padj"]
    .unstack("score_type")
    # .fillna(1)
    .sort_values(snakemake.wildcards["feature_set"])
    .merge(protein_coding_genes_df.rename(columns={"gene_id": "gene"}), on="gene", how="left")
    .set_index(["gene", "gene_name"])
    # # filter for Genebass genes
    # .dropna(subset=["Genebass"])
)
combined_regression_results_df

# %% [markdown]
# ## Plot

# %% [markdown]
# ## full plot

# %%
plot_df = combined_regression_results_df

x="Genebass"
y=snakemake.wildcards["feature_set"]

cutoff = 0.05
crop_pvalue = 10 ** -10

counts_x = (plot_df[x] < cutoff).sum()
counts_y = (plot_df[y] < cutoff).sum()

min_pval = min(
    plot_df[x].min(),
    plot_df[y].min(),
)

plot = (
    pn.ggplot(plot_df.fillna(1), pn.aes(x=x, y=y))
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
display(plot)

# %%
path = snakemake.params["output_basedir"] + "/pvalue_comp"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [markdown]
# ## cropped plot

# %%
plot_df = combined_regression_results_df

x="Genebass"
y=snakemake.wildcards["feature_set"]

cutoff = 0.05
crop_pvalue = 10 ** -10

plot_df = (
    plot_df
    .assign(**{
        x: np.fmax(crop_pvalue, plot_df[x].fillna(1)),
        y: np.fmax(crop_pvalue, plot_df[y].fillna(1)),
    })
    #.loc[(combined_regression_results_df[x] < cutoff) | (combined_regression_results_df[y] < cutoff), [x, y]]
)

counts_x = (plot_df[x] < cutoff).sum()
counts_y = (plot_df[y] < cutoff).sum()

min_pval = min(
    plot_df[x].min(),
    plot_df[y].min(),
)

plot = (
    pn.ggplot(plot_df, pn.aes(x=x, y=y))
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
display(plot)

# %%
path = snakemake.params["output_basedir"] + "/pvalue_comp_cropped"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [markdown]
# # save stats

# %%
significant_genes = (
    combined_regression_results_df
    .assign(**{
        f"{x}_signif": (combined_regression_results_df[x] < cutoff),
        f"{y}_signif": (combined_regression_results_df[y] < cutoff),
    })
    .query(f"{x}_signif | {y}_signif")
    .loc[:, [
        y,
        f"{y}_signif",
        x,
        f"{x}_signif",
    ]]
    .sort_values(y)
)

# %%
stats_df = (
    regression_results_df
    .assign(rsquared_diff=regression_results_df["rsquared"] - regression_results_df["rsquared_restricted"])
    .set_index("gene")
    .loc[:, [
        'n_observations', 'loglikelihood',
        'lr_stat', 'lr_df_diff', 'rsquared_restricted',
        'rsquared', 'rsquared_diff', 'lr_pval', 'padj',
    ]]
    .join(significant_genes[["Genebass", "Genebass_signif"]], how="right")
    .sort_values("padj")
)

# %%
with pd.option_context('display.float_format', '{:,.2g}'.format):
    display(
        stats_df
    )

# %%
stats_df.reset_index().to_parquet(snakemake.output["significant_genes_pq"], index=False)

# %%
stats_df.reset_index().to_csv(snakemake.output["significant_genes_tsv"], index=False, header=True, sep="\t", float_format='{:,.2g}'.format)

# %%

