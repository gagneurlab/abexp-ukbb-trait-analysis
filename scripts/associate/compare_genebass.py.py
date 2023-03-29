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

import polars as pl


# %%
import plotnine as pn
# import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt

# %%
from rep.notebook_init import setup_plot_style
setup_plot_style()

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

# %%
import matplotlib_venn

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
        rule_name = 'associate__compare_genebass',
        default_wildcards={
            # "phenotype_col": "Asthma",
            # "phenotype_col": "severe_LDL",
            # "phenotype_col": "triglycerides_f30870_0_0",
            #"phenotype_col": "hdl_cholesterol_f30760_0_0",
            "phenotype_col": "HDL_cholesterol",
            "feature_set": "LOFTEE_pLoF",
            # "feature_set": "AbExp_all_tissues",
            # "feature_set": "max_AbExp",
            # "covariates": "sex_age_genPC",
            "covariates": "sex_age_genPC_CLMP_PRS",
        }
    )

# %% {"tags": []}
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
if "plot_dpi" in snakemake.params:
    DPI = snakemake.params["plot_dpi"]
else:
    DPI=450

# %%
pval_cutoff = snakemake.params["pval_cutoff"]
pval_cutoff

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
protein_coding_genes_df = pl.scan_parquet(snakemake.input["protein_coding_genes_pq"])
protein_coding_genes_df.head().collect()

# %% [markdown] {"tags": []}
# ## read association results

# %%
snakemake.input["associations_pq"]

# %%
regression_results_df = (
    pl.scan_parquet(snakemake.input["associations_pq"] + "/*.parquet")
    .filter(pl.col("restricted_model_converged"))
    .filter(pl.col("full_model_converged"))
    .sort("rsquared", reverse=True)
    .drop([
        "term_pvals",
        "params",
    ])
    .with_column(pl.min([
        pl.col("lr_pval") * pl.count(),
        1.0,
    ]).alias("padj"))
    .with_column((pl.col("rsquared") - pl.col("rsquared_restricted")).alias("rsquared_diff"))
    .collect()
    # .to_pandas()
)

print(f"Corrected for {regression_results_df.shape[0]} association tests...")
# regression_results_df = regression_results_df.assign(padj=np.fmin(regression_results_df["lr_pval"] * regression_results_df.shape[0], 1))

display(regression_results_df)

# %%
regression_results_df.filter(pl.col("lr_stat") < -300)["gene"].to_pandas()

# %%
did_not_converge = (
    pl.scan_parquet(snakemake.input["associations_pq"] + "/*.parquet")
    .filter(
        (~ pl.col("restricted_model_converged")) | (~ pl.col("full_model_converged"))
    )
    .collect()
)
did_not_converge

# %%
did_not_converge["gene"].to_list()

# %%
regression_results_df["padj"].to_pandas().quantile(q=[0.5, 0.05, 0.001])

# %%
(regression_results_df["padj"] < 0.05).sum()

# %%
# regression_results_df[regression_results_df["padj"] < 0.05]

# %% {"tags": []}
genebass_300k_df = (
    pl.scan_parquet(snakemake.input["genebass_300k_pq"] + "/*.parquet")
    .filter(pl.col("annotation") == pl.lit("pLoF"))
    .filter(pl.col("phenocode") == pl.lit(phenocode))
    .rename({
        "gene_id": "gene",
        "Pvalue": "lr_pval",
    })
    .with_column(pl.min([
        pl.col("lr_pval") * pl.count(),
        1.0,
    ]).alias("padj"))
    .collect()
)
genebass_300k_df.schema

# %%
genebass_300k_df

# %% {"tags": []}
genebass_500k_df = (
    pl.scan_parquet(snakemake.input["genebass_500k_pq"] + "/*.parquet")
    .filter(pl.col("annotation") == pl.lit("pLoF"))
    .filter(pl.col("phenocode") == pl.lit(phenocode))
    .rename({
        "gene_id": "gene",
        "Pvalue": "lr_pval",
    })
    .with_column(pl.min([
        pl.col("lr_pval") * pl.count(),
        1.0,
    ]).alias("padj"))
    .collect()
)
genebass_500k_df.schema

# %%
genebass_500k_df

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

# %%
combined_regression_results_df = pd.concat([
    regression_results_df.with_column(pl.lit(snakemake.wildcards["feature_set"]).alias("score_type")).to_pandas(),
    genebass_300k_df.with_column(pl.lit("Genebass (300k WES)").alias("score_type")).to_pandas(),
    genebass_500k_df.with_column(pl.lit("Genebass (500k WES)").alias("score_type")).to_pandas(),
], join="inner")

combined_regression_results_df = (
    combined_regression_results_df
    .set_index(["gene", "score_type"])["padj"]
    .unstack("score_type")
    # .fillna(1)
    .sort_values(snakemake.wildcards["feature_set"])
    .merge(protein_coding_genes_df.collect().to_pandas().rename(columns={"gene_id": "gene"}), on="gene", how="left")
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
def plot_pvalues(plot_df, x, y, cutoff=0.05, crop_pvalue=10**-150):
    plot_df = (
        plot_df
        .assign(**{
            x: (np.fmax(crop_pvalue, plot_df[x].fillna(1))),
            y: (np.fmax(crop_pvalue, plot_df[y].fillna(1))),
        })
        #.loc[(combined_regression_results_df[x] < cutoff) | (combined_regression_results_df[y] < cutoff), [x, y]]
    )

    counts_x = (plot_df[x] < cutoff).sum()
    counts_y = (plot_df[y] < cutoff).sum()

    # min_pval = min(
    #     plot_df[x].min(),
    #     plot_df[y].min(),
    # )

    plot_df = (
        plot_df
        .assign(**{
            x: -np.log10(plot_df[x]),
            y: -np.log10(plot_df[y]),
        })
        #.loc[(combined_regression_results_df[x] < cutoff) | (combined_regression_results_df[y] < cutoff), [x, y]]
    )
    max_logpval = max(
        plot_df[x].max(),
        plot_df[y].max(),
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
        # + pn.facet_grid("~ Genebass_version")
        # + pn.scale_x_log10(limits=(min_pval, 1))
        # + pn.scale_y_log10(limits=(min_pval, 1))
        + pn.scale_x_continuous(limits=(0, max_logpval))
        + pn.scale_y_continuous(limits=(0, max_logpval))
        + pn.labs(
            title=f"Association between genes and '{phenotype_col}'\n(p-values, alpha={cutoff})",
            x=f"-log10('{x}' p-value)\n(n_signif={counts_x})",
            y=f"-log10('{y}' p-value)\n(n_signif={counts_y})",
        )
        + pn.theme(
            # axis_text_x=pn.element_text(rotation=90),
            # figure_size=(8, 8),
            title=pn.element_text(linespacing=1.4),
        )
        # + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
    )
    return plot_df, plot


# %%
plot_df, plot = plot_pvalues(
    combined_regression_results_df,
    x = "Genebass (300k WES)",
    y = snakemake.wildcards["feature_set"],
    cutoff = pval_cutoff,
    crop_pvalue = 10 ** -150,
)

# %%
plot

# %%
path = snakemake.params["output_basedir"] + "/pvalue_comp.Genebass_300k"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %%
plot_df, plot = plot_pvalues(
    combined_regression_results_df,
    x = "Genebass (500k WES)",
    y = snakemake.wildcards["feature_set"],
    cutoff = pval_cutoff,
    crop_pvalue = 10 ** -150,
)

# %%
plot

# %%
path = snakemake.params["output_basedir"] + "/pvalue_comp.Genebass_500k"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [markdown]
# ## cropped plot

# %%
plot_df, plot = plot_pvalues(
    combined_regression_results_df,
    x = "Genebass (300k WES)",
    y = snakemake.wildcards["feature_set"],
    cutoff = pval_cutoff,
    crop_pvalue = 10 ** -10,
)

# %%
plot

# %%
path = snakemake.params["output_basedir"] + "/pvalue_comp_cropped.Genebass_300k"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %%
plot_df, plot = plot_pvalues(
    combined_regression_results_df,
    x = "Genebass (500k WES)",
    y = snakemake.wildcards["feature_set"],
    cutoff = pval_cutoff,
    crop_pvalue = 10 ** -10,
)

# %%
plot

# %%
path = snakemake.params["output_basedir"] + "/pvalue_comp_cropped.Genebass_500k"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [markdown]
# ## Venn diagram

# %%
fig, ax = plt.subplots()
matplotlib_venn.venn3(
    (
        set(combined_regression_results_df[combined_regression_results_df[snakemake.wildcards["feature_set"]] < pval_cutoff].index.get_level_values("gene")),
        set(combined_regression_results_df[combined_regression_results_df["Genebass (300k WES)"] < pval_cutoff].index.get_level_values("gene")),
        set(combined_regression_results_df[combined_regression_results_df["Genebass (500k WES)"] < pval_cutoff].index.get_level_values("gene")),
    ),
    set_labels = (snakemake.wildcards["feature_set"], "Genebass (300k WES)", "Genebass (500k WES)"),
    ax=ax
)
display(ax)

# %%
snakemake.params["output_basedir"] + "/genebass_overlap"

# %%
path = snakemake.params["output_basedir"] + "/genebass_overlap"
fig.savefig(path + ".png", dpi=DPI)
fig.savefig(path + ".pdf", dpi=DPI)

# %% [raw]
# from IPython.display import Image
# display(Image(snakemake.params["output_basedir"] + "/genebass_overlap.png"))

# %% [markdown]
# # save stats

# %%
y = snakemake.wildcards["feature_set"]

significant_genes = (
    combined_regression_results_df
    .assign(**{
        "Genebass_300k_signif": (combined_regression_results_df["Genebass (300k WES)"] < pval_cutoff),
        "Genebass_500k_signif": (combined_regression_results_df["Genebass (500k WES)"] < pval_cutoff),
        f"{y}_signif": (combined_regression_results_df[y] < pval_cutoff),
    })
    .query(f"Genebass_300k_signif | Genebass_500k_signif | {y}_signif")
    # .query(f"{y}_signif")
    .loc[:, [
        y,
        f"{y}_signif",
        "Genebass (300k WES)",
        "Genebass_300k_signif",
        "Genebass (500k WES)",
        "Genebass_500k_signif",
    ]]
    .sort_values(y)
)

# %%
stats_df = (
    regression_results_df.to_pandas()
    .set_index("gene")
    .loc[:, [
        'n_observations',
        'restricted_model_converged',
        'full_model_converged',
        'restricted_model_llf',
        'full_model_llf',
        'lr_stat', 'lr_df_diff',
        'rsquared_restricted',
        "rsquared_restricted_raw",
        'rsquared',
        "rsquared_raw",
        'rsquared_diff',
        'lr_pval',
        'padj',
    ]]
    .join(significant_genes[[
        "Genebass (300k WES)",
        "Genebass_300k_signif",
        "Genebass (500k WES)",
        "Genebass_500k_signif",
    ]], how="right")
    .sort_values("padj")
    .assign(**snakemake.wildcards)
)

# %%
with pd.option_context('display.float_format', '{:,.2g}'.format):
    display(
        stats_df
    )

# %%
stats_df.reset_index().to_parquet(snakemake.output["significant_genes_pq"], index=False)

# %%
stats_df.reset_index().to_csv(snakemake.output["significant_genes_tsv"], index=False, header=True, sep="\t")

# %%

