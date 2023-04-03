# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
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
        rule_name = 'associate__compare_monti',
        default_wildcards={
            # "phenotype_col": "Asthma",
            # "phenotype_col": "severe_LDL",
            # "phenotype_col": "triglycerides_f30870_0_0",
            # "phenotype_col": "Alkaline_phosphatase",
            "phenotype_col": "HDL_cholesterol",
            # "feature_set": "LOFTEE_pLoF",
            "feature_set": "AbExp_all_tissues",
            # "feature_set": "max_AbExp",
            # "covariates": "sex_age_genPC",
            "covariates": "sex_age_genPC_CLMP_PRS",
        }
    )

# %%
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

# %% [markdown]
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

# %% [markdown]
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
    .with_columns([
        pl.min([
            pl.col("lr_pval") * pl.count(),
            1.0,
        ]).alias("padj"),
        (pl.col("rsquared") - pl.col("rsquared_restricted")).alias("rsquared_diff"),
    ])
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

# %% [markdown]
# ### genebass

# %%
genebass_300k_df = (
    pl.scan_parquet(snakemake.input["genebass_300k_pq"] + "/*.parquet")
    .filter(pl.col("annotation") == pl.lit("pLoF"))
    .filter(pl.col("phenocode") == pl.lit(phenocode))
    .rename({
        "gene_id": "gene",
        "Pvalue": "lr_pval",
    })
    .with_columns(pl.min([
        pl.col("lr_pval") * pl.count(),
        1.0,
    ]).alias("padj"))
    .collect()
)
genebass_300k_df.schema

# %%
genebass_300k_df

# %%
genebass_500k_df = (
    pl.scan_parquet(snakemake.input["genebass_500k_pq"] + "/*.parquet")
    .filter(pl.col("annotation") == pl.lit("pLoF"))
    .filter(pl.col("phenocode") == pl.lit(phenocode))
    .rename({
        "gene_id": "gene",
        "Pvalue": "lr_pval",
    })
    .with_columns(pl.min([
        pl.col("lr_pval") * pl.count(),
        1.0,
    ]).alias("padj"))
    .collect()
)
genebass_500k_df.schema

# %%
genebass_500k_df

# %%
genebass_significant_genes = (
    genebass_500k_df
    .select(
        "gene",
        (pl.col("padj") < pl.lit(0.05)).alias("genebass_500k_significant")
    )
)
genebass_significant_genes

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

# %% [markdown]
# ### monti

# %%
# pl.scan_parquet(snakemake.input["monti_results_pq"]).select("pheno").unique().sort("pheno").collect()["pheno"].to_list()

# %%
monti_mapping_df = pl.read_csv(snakemake.input["monti_mapping_csv"])
monti_mapping_df

monti_results_df = (
    pl.scan_parquet(snakemake.input["monti_results_pq"])
    .rename({"pheno": "monti_pheno"})
)
monti_results_df = (
    monti_mapping_df.lazy()
    .join(monti_results_df, on="monti_pheno", how="left")
    .drop("monti_pheno")
    .filter(pl.col("phenotype") == pl.lit(phenotype_col))
    .select([
        'phenotype',
        pl.col('ensembl97_id').alias("gene"),
        'gene_name',
        # 'N',
        pl.min([
            "pv_score_linb_pLOF",
            "pv_slrt_linwb_cct_miss",
            "pv_slrt_linwcollapsed_cct_miss"
        ]).alias("lr_pval"),
        #'pv_score_linb_pLOF',
        # 'nCarrier_pLOF',
        # 'n_snp_pLOF',
        # 'nCarrier_miss',
        # 'cumMAC_miss',
        # 'n_snp_miss',
        # 'n_cluster_miss',
        # 'pv_slrt_linwb_cct_miss',
        # 'pv_slrt_linwcollapsed_cct_miss',
        # 'nCarrier_splice',
        # 'nCarrier_notLOF_splice',
        # 'cumMAC_splice',
        # 'n_snp_splice',
        # 'pv_slrt_linwb_cct_splice',
        # 'pv_slrt_linw_cct_splice',
        # 'nCarrier_rbp',
        # 'cumMAC_rbp',
        # 'n_snp_rbp',
        # 'n_snp_notLOF_rbp',
        # 'pv_lrt_linwcholesky_rbp'
    ])
    .filter(pl.col("lr_pval").is_not_null())
    .sort("lr_pval")
    .with_columns([
        pl.min([
            pl.col("lr_pval") * pl.count(),
            1.0,
        ]).alias("padj")
    ])
    .collect()
)
monti_results_df

# %% [markdown]
# # join everything

# %%
target_cols = ["phenotype", "gene", "score_type", "lr_pval", "padj"]

combined_regression_results_df = pl.concat([
    (
        regression_results_df
        .with_columns([
            pl.lit(snakemake.wildcards["phenotype_col"]).alias("phenotype"),
            pl.lit(snakemake.wildcards["feature_set"]).alias("score_type"),
        ])
        .select(target_cols)
    ),
    #(
    #    genebass_300k_df
    #    .with_columns([
    #        pl.lit(snakemake.wildcards["phenotype_col"]).alias("phenotype"),
    #        pl.lit("Genebass (300k WES)").alias("score_type"),
    #    ])
    #    .select(target_cols)
    #),
    #(
    #    genebass_500k_df
    #    .with_columns([
    #        pl.lit(snakemake.wildcards["phenotype_col"]).alias("phenotype"),
    #        pl.lit("Genebass (500k WES)").alias("score_type"),
    #    ])
    #    .select(target_cols)
    #),
    (
        monti_results_df
        .with_columns(pl.lit("Monti").alias("score_type"))
        .select(target_cols)
    ),
])
combined_regression_results_df = (
    genebass_significant_genes
    .join(
        combined_regression_results_df.select("score_type").unique(),
        how="cross"
    )
    .join(
        combined_regression_results_df.select("phenotype").unique(),
        how="cross"
    )
    .join(
        combined_regression_results_df,
        on=["gene", "phenotype", "score_type"],
        how="left"
    )
    .with_columns([
        pl.col("lr_pval").fill_null(1).fill_nan(1),
        pl.col("padj").fill_null(1).fill_nan(1),
    ])
    .sort("score_type", "padj")
)
# compute ranks
combined_regression_results_df = (
    combined_regression_results_df
    .groupby("phenotype", "score_type")
    .apply(lambda group: group.with_columns(pl.col("lr_pval").rank().alias("rank")))
)
## filter for monti genes
#combined_regression_results_df = (
#    combined_regression_results_df
#    .filter(pl.col("gene").is_in(monti_results_df["gene"]))
#)
combined_regression_results_df

# %% [markdown]
# ## Plot

# %%
plot_df = (
    combined_regression_results_df
    .groupby("phenotype", "score_type")
    .apply(lambda group: (
        group
        .sort("rank")
        .with_columns(pl.col("genebass_500k_significant").cumsum().alias("n_true"))
    ))
    .with_columns([
        pl.lit(v).alias(k) for k, v in snakemake.wildcards.items()
    ])
)
plot_df

# %%
plot = (
    pn.ggplot(plot_df, pn.aes(x="rank", y="n_true", color="score_type"))
    + pn.geom_line()
)
plot

# %%
path = snakemake.params["output_basedir"] + "/compare_monti"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [raw]
# from IPython.display import Image
# display(Image(snakemake.params["output_basedir"] + "/genebass_overlap.png"))

# %%
plot_df.write_parquet(snakemake.output["compare_monti_pq"], compression="snappy", statistics=True, use_pyarrow=True)

# %%

