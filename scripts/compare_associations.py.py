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
from rep.notebook_init import setup_plot_style
setup_plot_style()

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
        rule_name = 'compare_associations',
        default_wildcards={
            "comparison": "all",
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
# # Load configuration

# %%
with open(snakemake.params["config_yaml"], "r") as fd:
    config = yaml.safe_load(fd)

# %%
print(json.dumps(config, indent=2, default=str))

# %%
pval_cutoff = config["pval_cutoff"]
pval_cutoff

# %% [markdown] {"tags": []}
# # Read features

# %% [markdown] {"tags": []}
# ## read association results

# %%
snakemake.input["associations_pq"]

# %%
regression_results_df = (
    spark.read.parquet(*snakemake.input["associations_pq"])
    # .sort("rsquared", reverse=True)
    # .drop([
    #     "term_pvals",
    #     "params",
    # ])
    # .with_column(pl.min([
    #     pl.col("lr_pval") * pl.count(),
    #     1.0,
    # ]).alias("padj"))
    # .with_column((pl.col("rsquared") - pl.col("rsquared_restricted")).alias("rsquared_diff"))
    # .collect()
    # # .to_pandas()
)

# print(f"Corrected for {regression_results_df.shape[0]} association tests...")
# regression_results_df = regression_results_df.assign(padj=np.fmin(regression_results_df["lr_pval"] * regression_results_df.shape[0], 1))

regression_results_df.printSchema()

# %%
grouping = ["phenotype_col", "covariates", "feature_set"]

# %%
counts = (
    regression_results_df
    .groupby(grouping)
    .count()
    .withColumnRenamed("count", "num_association_tests")
    .sort(grouping)
    .persist()
)
counts.toPandas()

# %%
adj_regression_results_df = (
    regression_results_df
    .join(counts, on=grouping, how="left")
    .withColumn("padj", f.col("lr_pval") * f.col("num_association_tests"))
    .withColumn("is_significant", f.col("padj") < f.lit(pval_cutoff))
    .withColumn("is_significant", f.when(f.col("is_significant").isNotNull(), f.col("is_significant")).otherwise(f.lit(False)))
    .withColumn("rsquared_diff", f.col("rsquared") - f.col("rsquared_restricted"))
    .sort([
        *grouping,
        "rsquared_diff",
    ], ascending=False)
)
adj_regression_results_df.printSchema()

# %%
significant_associations = (
    adj_regression_results_df
    .filter(f.col("padj") < f.lit(pval_cutoff))
    .sort([
        *grouping,
        "rsquared_diff",
    ], ascending=False)
    .persist()
)

# %%
(
    significant_associations
    .write
    .parquet(snakemake.output["significant_genes_pq"], mode="overwrite")
)

# %%
adj_regression_results_pd_df = (
    adj_regression_results_df
    .drop("term_pvals", "params")
    .toPandas()
)
adj_regression_results_pd_df

# %%
significant_associations_pd_df = (
    significant_associations
    .drop("term_pvals", "params")
    .toPandas()
)
significant_associations_pd_df

# %%
with pd.option_context('display.max_rows', 500, 'display.max_columns', 20):
    display(
        significant_associations_pd_df
        .sort_values(["phenotype_col", "gene", "padj"])
        .drop(columns=[
            "n_observations",
            "num_association_tests",
        ])
        .set_index(["gene", "phenotype_col"])
    )

# %%
num_significant_associations = (
    significant_associations
    .groupby(grouping)
    .count()
    .toPandas()
)
num_significant_associations

# %%
num_significant_associations = (
    adj_regression_results_df
    .groupby(grouping)
    .agg(
        f.sum(f.col("is_significant").cast(t.IntegerType())).alias("num_significant"),
        f.count(f.col("is_significant")).alias("total_association_tests"),
    )
    .toPandas()
)
num_significant_associations

# %%
path = snakemake.params["output_basedir"] + "/num_significant_associations"
num_significant_associations.to_csv(f"{path}.csv", index=False)
num_significant_associations.to_parquet(f"{path}.parquet", index=False)

# %% [markdown]
# ## Plot

# %% [markdown]
# ## boxplot

# %%
plot_df = adj_regression_results_pd_df.fillna({
    "padj": 1.,
    "lr_pval": 1.,
    "rsquared_diff": 0.,
})

plot_df = plot_df.assign(covariates=plot_df["covariates"].str.replace("_", "\n+ "))

# crop_pvalue = 10 ** -10

plot = (
    pn.ggplot(plot_df, pn.aes(x="feature_set", y="rsquared_diff"))
    # + pn.geom_abline(slope=1, color="black", linetype="dashed")
    # + pn.geom_hline(yintercept=0.05, color="red", linetype="dashed")
    # + pn.geom_vline(xintercept=0.05, color="red", linetype="dashed")
    # + pn.geom_rect(
    #     xmax=1,
    #     xmin=0.05,
    #     ymax=1,
    #     ymin=0.05,
    #     color="red", 
    #     linetype="dashed"
    # )
    # + pn.geom_point(data=plot_df.query("is_significant"), color="red")
    + pn.geom_boxplot()
    + pn.geom_point(pn.aes(color="is_significant"), data=plot_df.query("is_significant"))
    # + pn.geom_text(pn.aes(label="label", y=-0.1), data=num_significant_associations.assign(
    #     label=num_significant_associations["num_significant"].apply(lambda s: f"n_signif={s}"),
    #     covariates=num_significant_associations["covariates"].str.replace("_", "\n+ ")
    # ))
    # + pn.geom_point()
    #+ pn.coord_fixed()
    # + pn.scale_x_log10(limits=(min_pval, 1))
    # + pn.scale_y_log10(limits=(min_pval, 1))
    + pn.labs(
        title=f"Difference in explained variance (r²) for gene-trait associations\nwhen adding a certain feature set",
        # f"\n(p-values, alpha={cutoff})",
        x=f"feature set",
        y=f"r²(full) - r²(restricted)",
    )
    # + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
    + pn.theme(axis_text_x=pn.element_text(
        rotation=90,
        # hjust=1
    ))
    + pn.facet_grid("phenotype_col ~ covariates")
    # + pn.coord_flip()
)

# %%
display(plot)

# %%
snakemake.params["output_basedir"]

# %%
path = snakemake.params["output_basedir"] + "/rsquared_diff"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [markdown]
# ## scatter plot r²-diff

# %%
import itertools

keys = adj_regression_results_pd_df["feature_set"].unique().tolist()

plot_df = adj_regression_results_pd_df.set_index([*grouping, "gene"]).unstack("feature_set").fillna({
    "is_significant": False,
    "padj": 1.,
    "lr_pval": 1.,
    "rsquared_diff": 0.,
})
plot_df.columns = [f"{b}.{a}" for a, b in plot_df.columns]

for feature_x, feature_y in list(itertools.combinations(keys, 2)):
    print(f"Plotting '{feature_x}' vs '{feature_y}'...")
    
    signif_map = {
        (True, True): "both",
        (True, False): feature_x,
        (False, True): feature_y,
        (False, False): "none",
    }
    
    sub_plot_df = (
        plot_df
        .reset_index()
        .assign(
            is_significant=[
                signif_map[(v[0], v[1])] for k, v in plot_df[[
                    f'{feature_x}.is_significant',
                    f'{feature_y}.is_significant'
                ]].iterrows()
            ]
        )
    )

    plot = (
        pn.ggplot(sub_plot_df, pn.aes(
            x=f'{feature_x}.rsquared_diff',
            y=f'{feature_y}.rsquared_diff',
            color="is_significant",
        ))
        + pn.geom_point()
        + pn.ggtitle(f"Comparison between '{feature_x}' and '{feature_y}'\n(r²(full) - r²(restricted))")
        + pn.labs(
            x=f"{feature_x}\n(r²(full) - r²(restricted))",
            y=f"{feature_y}\n(r²(full) - r²(restricted))",
        )
        + pn.facet_grid("phenotype_col ~ covariates")
        + pn.theme(axis_text_x=pn.element_text(
            rotation=30,
            # hjust=1
        ))
    )

    
    path = snakemake.params["output_basedir"] + f"/feature_comp.r2_diff@{feature_x}__vs__{feature_y}"
    
    print(f"Saving to '{path}'")
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)



# %% [markdown]
# ## scatter plot p-values

# %%
def scatter_plot_pval(plot_df, x, y, label_x, label_y, cutoff=pval_cutoff, crop_pvalue=None):
    if crop_pvalue is not None:
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
        + pn.scale_x_log10()
        + pn.scale_y_log10()
        #+ pn.coord_fixed()
        + pn.labs(
            x=f"{label_x}\n(n_signif={counts_x})",
            y=f"{label_y}\n(n_signif={counts_y})",
        )
        # + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
    )
    
    if crop_pvalue is not None:
        plot = (
            plot
            + pn.scale_x_log10(limits=(min_pval, 1))
            + pn.scale_y_log10(limits=(min_pval, 1))
        )
    
    return plot


# %%
import itertools

keys = adj_regression_results_pd_df["feature_set"].unique().tolist()

plot_df = adj_regression_results_pd_df.set_index([*grouping, "gene"]).unstack("feature_set").fillna({
    "is_significant": False,
    "padj": 1.,
    "lr_pval": 1.,
    "rsquared_diff": 0.,
})
plot_df.columns = [f"{b}.{a}" for a, b in plot_df.columns]

for feature_x, feature_y in list(itertools.combinations(keys, 2)):
    print(f"Plotting '{feature_x}' vs '{feature_y}'...")
    
    signif_map = {
        (True, True): "both",
        (True, False): feature_x,
        (False, True): feature_y,
        (False, False): "none",
    }
    
    sub_plot_df = (
        plot_df
        .reset_index()
        .assign(
            is_significant=[
                signif_map[(v[0], v[1])] for k, v in plot_df[[
                    f'{feature_x}.is_significant',
                    f'{feature_y}.is_significant'
                ]].iterrows()
            ]
        )
    )
    
    plot = (
        scatter_plot_pval(
            sub_plot_df,
            x=f'{feature_x}.padj',
            y=f'{feature_y}.padj',
            label_x=feature_x,
            label_y=feature_y,
        )
        + pn.ggtitle(f"Comparison between '{feature_x}' and '{feature_y}'\n(p-values, alpha={pval_cutoff})")
        + pn.facet_grid("phenotype_col ~ covariates")
        + pn.theme(axis_text_x=pn.element_text(
            rotation=30,
            # hjust=1
        ))
    )
    
    path = snakemake.params["output_basedir"] + f"/feature_comp.padj@{feature_x}__vs__{feature_y}"
    
    print(f"Saving to '{path}'")
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)



# %% [markdown]
# ## scatter plot p-values (cropped)

# %%
import itertools

keys = adj_regression_results_pd_df["feature_set"].unique().tolist()

plot_df = adj_regression_results_pd_df.set_index([*grouping, "gene"]).unstack("feature_set").fillna({
    "is_significant": False,
    "padj": 1.,
    "lr_pval": 1.,
    "rsquared_diff": 0.,
})
plot_df.columns = [f"{b}.{a}" for a, b in plot_df.columns]

for feature_x, feature_y in list(itertools.combinations(keys, 2)):
    print(f"Plotting '{feature_x}' vs '{feature_y}'...")
    
    signif_map = {
        (True, True): "both",
        (True, False): feature_x,
        (False, True): feature_y,
        (False, False): "none",
    }
    
    sub_plot_df = (
        plot_df
        .reset_index()
        .assign(
            is_significant=[
                signif_map[(v[0], v[1])] for k, v in plot_df[[
                    f'{feature_x}.is_significant',
                    f'{feature_y}.is_significant'
                ]].iterrows()
            ]
        )
    )
    
    plot = (
        scatter_plot_pval(
            sub_plot_df,
            x=f'{feature_x}.padj',
            y=f'{feature_y}.padj',
            label_x=feature_x,
            label_y=feature_y,
            crop_pvalue=10 ** -10,
        )
        + pn.ggtitle(f"Comparison between '{feature_x}' and '{feature_y}'\n(p-values, alpha={pval_cutoff})")
        + pn.facet_grid("phenotype_col ~ covariates")
        + pn.theme(axis_text_x=pn.element_text(
            rotation=30,
            # hjust=1
        ))
    )
    
    path = snakemake.params["output_basedir"] + f"/feature_comp.padj_cropped@{feature_x}__vs__{feature_y}"
    
    print(f"Saving to '{path}'")
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)



# %%
