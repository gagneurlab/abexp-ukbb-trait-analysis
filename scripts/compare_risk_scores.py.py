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

# import glow


# %%
import plotnine as pn
# import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt

# %%
import textwrap

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
spark = init_spark(enable_glow=False)

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
        rule_name = 'compare_risk_scores',
        default_wildcards={
            "comparison": "all", 
            # "comparison": "paper_figure", 
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
# ## read PRS results

# %% {"tags": []}
snakemake.input["predictions_pq"]

# %%
predictions_df = (
    spark.read.parquet(*snakemake.input["predictions_pq"])
)
predictions_df.printSchema()

# %% [markdown]
# ## compute r2 scores

# %%
from sklearn.metrics import r2_score

# %%
r2_score_df = []
for path in snakemake.input["predictions_pq"]:
    df = pd.read_parquet(path)
    r2_score_df.append(
        pd.DataFrame({
            "rsquared_full": [r2_score(y_true=df["measurement"], y_pred=df["full_model_pred"])],
            "rsquared_restricted": [r2_score(y_true=df["measurement"], y_pred=df["restricted_model_pred"])],
            "rsquared_basic": [r2_score(y_true=df["measurement"], y_pred=df["basic_model_pred"])],
            **{c: [df[c].iloc[0]] for c in [
                "phenotype_col",
                "feature_set",
                "covariates",
            ]},
        })
    )
r2_score_df = pd.concat(r2_score_df).reset_index(drop=True)
r2_score_df

# %%
path = snakemake.params["output_basedir"] + f"/rsquared"
r2_score_df.to_csv(path + ".csv", index=False)
r2_score_df.to_parquet(path + ".parquet", index=False)


# %% [markdown]
# ## compute ranks

# %%
def assign_ranks(df: pd.DataFrame):
    df = df.assign(**{
        "measurement_rank": df["measurement"].rank(pct=True),
        "full_model_pred_rank": df["full_model_pred"].rank(pct=True),
        "restricted_model_pred_rank": df["restricted_model_pred"].rank(pct=True),
        "basic_model_pred_rank": df["basic_model_pred"].rank(pct=True),
        "total": df.shape[0],
    }).astype({
        "measurement_rank": "float64",
        "full_model_pred_rank": "float64",
        "restricted_model_pred_rank": "float64",
        "basic_model_pred_rank": "float64",
        "total": "int64",
    })
    return df


# %%
from copy import deepcopy

# %%
returned_schema = (
    deepcopy(predictions_df.schema)
    .add("measurement_rank", t.DoubleType())
    .add("full_model_pred_rank", t.DoubleType())
    .add("restricted_model_pred_rank", t.DoubleType())
    .add("basic_model_pred_rank", t.DoubleType())
    .add("total", t.LongType())
)

# %%
transformed_predictions_df = (
    predictions_df
    .groupby(["phenotype_col", "feature_set", "covariates"])
    .applyInPandas(assign_ranks, schema=returned_schema)
)

# %%
transformed_predictions_df.printSchema()

# %%
target_quantiles = [0.01, 0.05, 0.1, 0.2]
# lower_upper = ["lower", "upper"]

target_quantiles_df = (
    pd.DataFrame({"measurement_quantile": target_quantiles})
    .merge(pd.DataFrame({"prediction_quantile": target_quantiles}), how='cross')
    # .merge(pd.DataFrame({"bound": lower_upper}), how='cross')
)
target_quantiles_df

# %%
target_quantiles_sdf = spark.createDataFrame(target_quantiles_df)

# %%
quantiles_df = (
    transformed_predictions_df.crossJoin(target_quantiles_sdf)
    .groupby([*target_quantiles_sdf.columns, "phenotype_col", "feature_set", "covariates"])
    .agg(
        f.struct([
            f.sum((f.col("measurement_rank") > (1 - f.col("measurement_quantile"))).cast(t.LongType())).alias("total"),
            f.sum((
                (f.col("full_model_pred_rank") > (1 - f.col("prediction_quantile"))) & (f.col("measurement_rank") > (1 - f.col("measurement_quantile")))
            ).cast(t.LongType())).alias("full_model"),
            f.sum((
                (f.col("restricted_model_pred_rank") > (1 - f.col("prediction_quantile"))) & (f.col("measurement_rank") > (1 - f.col("measurement_quantile")))
            ).cast(t.LongType())).alias("restricted_model"),
            f.sum((
                (f.col("basic_model_pred_rank") > (1 - f.col("prediction_quantile"))) & (f.col("measurement_rank") > (1 - f.col("measurement_quantile")))
            ).cast(t.LongType())).alias("basic_model"),
            f.lit("upper").alias("bound"),
        ]).alias("upper"),
        f.struct([
            f.sum((f.col("measurement_rank") < f.col("measurement_quantile")).cast(t.LongType())).alias("total"),
            f.sum((
                (f.col("full_model_pred_rank") < f.col("prediction_quantile")) & (f.col("measurement_rank") < f.col("measurement_quantile"))
            ).cast(t.LongType())).alias("full_model"),
            f.sum((
                (f.col("restricted_model_pred_rank") < f.col("prediction_quantile")) & (f.col("measurement_rank") < f.col("measurement_quantile"))
            ).cast(t.LongType())).alias("restricted_model"),
            f.sum((
                (f.col("basic_model_pred_rank") < f.col("prediction_quantile")) & (f.col("measurement_rank") < f.col("measurement_quantile"))
            ).cast(t.LongType())).alias("basic_model"),
            f.lit("lower").alias("bound"),
        ]).alias("lower")
    )
    .withColumn("counts", f.explode(f.array(f.col("upper"), f.col("lower"))))
    .select("*", "counts.*")
    .drop("upper", "lower", "counts")
)
quantiles_df.printSchema()

# %%
quantiles_pd_df = quantiles_df.toPandas()
quantiles_pd_df

# %%
path = snakemake.params["output_basedir"] + f"/quantile_counts"
quantiles_pd_df.to_csv(path + ".csv", index=False)
quantiles_pd_df.to_parquet(path + ".parquet", index=False)

# %% [markdown]
# ## Plot

# %% [markdown]
# ## scatter-plot r2

# %%
r2_score_df

# %%
plot_df = r2_score_df
grouping = ['phenotype_col', 'feature_set', 'covariates']
keys = plot_df["feature_set"].unique().tolist()

unstacked_plot_full_df = plot_df.set_index(grouping)["rsquared_full"].unstack("feature_set")
unstacked_plot_full_df
unstacked_plot_restricted_df = plot_df.set_index(grouping)["rsquared_restricted"].unstack("feature_set")
unstacked_plot_restricted_df

# %% [markdown]
# ### bar plot difference

# %%
import itertools

# list(itertools.combinations(keys, 2))
for feature_x, feature_y in list(itertools.product(keys, keys)):
    if feature_x == feature_y:
        continue
    
    subset_plot_df = unstacked_plot_full_df.reset_index()
    assert np.all(unstacked_plot_restricted_df[feature_x].values == unstacked_plot_restricted_df[feature_y].values)
    subset_plot_df = subset_plot_df.assign(**{
        f"difference_to_{feature_x}": subset_plot_df[feature_y] - subset_plot_df[feature_x],
        f"proportional_difference_to_{feature_x}": (subset_plot_df[feature_y] - subset_plot_df[feature_x]) / subset_plot_df[feature_x],
        "phenotype_col": subset_plot_df["phenotype_col"].str.replace("_", " "),
        "baseline": unstacked_plot_restricted_df[feature_x].values,
        # f"difference_to_baseline": subset_plot_df[feature_y] - unstacked_plot_restricted_df[feature_x].values,
    })
    subset_plot_df

    plot = (
        pn.ggplot(subset_plot_df.reset_index(), pn.aes(
            x=f"reorder(phenotype_col, difference_to_{feature_x})",
            y=f"difference_to_{feature_x}",
        ))
        + pn.geom_bar(
            stat="identity",
        )
        + pn.labs(
            x=f"""phenotype""",
            y=f"""difference in r² between '{feature_y.replace("_", " ")}' and '{feature_x.replace("_", " ")}'""",
            title=f"Comparison of phenotype prediction models using different feature sets",
        )
        + pn.theme(
            # legend_text=pn.element_text(linespacing=1.4),
            figure_size=(8, 12),
            axis_text_x=pn.element_text(
            #     rotation=45,
            #     hjust=1
                # vjust=10,
            ),
            # strip_text_y=pn.element_text(
            #     rotation=0,
            # ),
            title=pn.element_text(linespacing=1.4, vjust=-10),
            axis_title_x=pn.element_text(linespacing=1.4, vjust=-10),
        )
        # + pn.coord_equal()
        + pn.coord_flip()
    )
    display(plot)
    
    path = snakemake.params["output_basedir"] + f"/r2_bar_plot_difference.{feature_x}__vs__{feature_y}"
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)


# %% [markdown]
# ### bar plot proportional difference

# %%
import itertools
import mizani

# list(itertools.combinations(keys, 2))
for feature_x, feature_y in list(itertools.product(keys, keys)):
    if feature_x == feature_y:
        continue
    
    subset_plot_df = unstacked_plot_full_df.reset_index()
    assert np.all(unstacked_plot_restricted_df[feature_x].values == unstacked_plot_restricted_df[feature_y].values)
    subset_plot_df = subset_plot_df.assign(**{
        f"difference_to_{feature_x}": subset_plot_df[feature_y] - subset_plot_df[feature_x],
        f"proportional_difference_to_{feature_x}": (subset_plot_df[feature_y] - subset_plot_df[feature_x]) / subset_plot_df[feature_x],
        "phenotype_col": subset_plot_df["phenotype_col"].str.replace("_", " "),
        "baseline": unstacked_plot_restricted_df[feature_x].values,
        # f"difference_to_baseline": subset_plot_df[feature_y] - unstacked_plot_restricted_df[feature_x].values,
    })
    subset_plot_df

    plot = (
        pn.ggplot(subset_plot_df.reset_index(), pn.aes(
            x=f"reorder(phenotype_col, proportional_difference_to_{feature_x})",
            y=f"proportional_difference_to_{feature_x}",
        ))
        + pn.geom_bar(
            stat="identity",
        )
        + pn.scale_y_continuous(
            labels=mizani.formatters.percent_format()
        )
        + pn.labs(
            x=f"""phenotype""",
            y=f"""proportional difference in r² between '{feature_y.replace("_", " ")}' and '{feature_x.replace("_", " ")}'""",
            title=f"Comparison of phenotype prediction models using different feature sets",
        )
        + pn.theme(
            # legend_text=pn.element_text(linespacing=1.4),
            figure_size=(8, 12),
            axis_text_x=pn.element_text(
            #     rotation=45,
            #     hjust=1
                # vjust=10,
            ),
            # strip_text_y=pn.element_text(
            #     rotation=0,
            # ),
            title=pn.element_text(linespacing=1.4, vjust=-10),
            axis_title_x=pn.element_text(linespacing=1.4, vjust=-10),
        )
        # + pn.coord_equal()
        + pn.coord_flip()
    )
    display(plot)
    
    path = snakemake.params["output_basedir"] + f"/r2_bar_plot_proportional_difference.{feature_x}__vs__{feature_y}"
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)


# %% [markdown]
# ### scatter plot

# %%
import itertools

# list(itertools.combinations(keys, 2))
for feature_x, feature_y in list(itertools.product(keys, keys)):
    if feature_x == feature_y:
        continue
    
    subset_plot_df = unstacked_plot_full_df.reset_index()
    assert np.all(unstacked_plot_restricted_df[feature_x].values == unstacked_plot_restricted_df[feature_y].values)
    subset_plot_df = subset_plot_df.assign(**{
        f"difference_to_{feature_x}": subset_plot_df[feature_y] - subset_plot_df[feature_x],
        f"proportional_difference_to_{feature_x}": (subset_plot_df[feature_y] - subset_plot_df[feature_x]) / subset_plot_df[feature_x],
        "phenotype_col": subset_plot_df["phenotype_col"].str.replace("_", " "),
        "baseline": unstacked_plot_restricted_df[feature_x].values,
        # f"difference_to_baseline": subset_plot_df[feature_y] - unstacked_plot_restricted_df[feature_x].values,
    })
    subset_plot_df

    plot = (
        pn.ggplot(subset_plot_df.reset_index(), pn.aes(
            x=f"{feature_x}",
            y=f"{feature_y}",
            fill=f"phenotype_col",
        ))
        + pn.geom_point(
            # pn.aes(width=.95, height=.95),
            size=3,
        )
        + pn.geom_abline(slope=1, linetype="dashed")
        # + pn.scale_fill_gradient2(
        #     # low = muted("red"),
        #     # mid = "white",
        #     # high = muted("blue"),
        #     midpoint = 0,
        #     labels=mizani.formatters.percent_format()
        # )
        + pn.labs(
            x=f"""{feature_x.replace("_", " ")} (r²)""",
            y=f"""{feature_y.replace("_", " ")} (r²)""",
            title=f"Comparison of phenotype prediction models using different feature sets",
        )
        + pn.theme(
            legend_text=pn.element_text(linespacing=1.4),
            # figure_size=(8, 8),
            # axis_text_x=pn.element_text(
            #     rotation=45,
            #     hjust=1
            # ),
            # strip_text_y=pn.element_text(
            #     rotation=0,
            # ),
            title=pn.element_text(linespacing=1.4, vjust=-10),
        )
        + pn.coord_equal()
        # + pn.coord_flip()
    )
    display(plot)
    
    path = snakemake.params["output_basedir"] + f"/r2_scatter.{feature_x}__vs__{feature_y}"
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)


# %% [markdown]
# ### bar plot absolute r²

# %%
import itertools

# list(itertools.combinations(keys, 2))
for feature_x, feature_y in list(itertools.product(keys, keys)):
    if feature_x == feature_y:
        continue
    
    subset_plot_df = unstacked_plot_full_df.reset_index()
    assert np.all(unstacked_plot_restricted_df[feature_x].values == unstacked_plot_restricted_df[feature_y].values)
    subset_plot_df = subset_plot_df.assign(**{
        f"difference_to_{feature_x}": subset_plot_df[feature_y] - subset_plot_df[feature_x],
        "phenotype_col": subset_plot_df["phenotype_col"].str.replace("_", " "),
        "baseline": unstacked_plot_restricted_df[feature_x].values,
    })
    subset_plot_df = (
        subset_plot_df
        .reset_index()
        .melt(id_vars=["phenotype_col", "covariates"], value_vars=["baseline", feature_x, feature_y], value_name="rsquared")
    )
    subset_plot_df = (
        subset_plot_df
        .assign(**{
            "feature_set": subset_plot_df["feature_set"].str.replace("_", " "),
        })
    )
    subset_plot_df

    plot = (
        pn.ggplot(subset_plot_df, pn.aes(
            x=f"reorder(phenotype_col, rsquared)",
            y=f"rsquared",
            fill=f"feature_set",
        ))
        + pn.geom_bar(stat="identity", position="dodge")
        # + pn.scale_fill_gradient2(
        #     # low = muted("red"),
        #     # mid = "white",
        #     # high = muted("blue"),
        #     midpoint = 0,
        #     labels=mizani.formatters.percent_format()
        # )
        + pn.labs(
            x=f"""phenotype""",
            y=f"""r²""",
            title=f"Comparison of PRS models using different feature sets",
        )
        + pn.theme(
            legend_text=pn.element_text(linespacing=1.4),
            figure_size=(8, 12),
            # axis_text_x=pn.element_text(
            #     rotation=45,
            #     hjust=1
            # ),
            # strip_text_y=pn.element_text(
            #     rotation=0,
            # ),
            title=pn.element_text(linespacing=1.4, vjust=-10),
        )
        # + pn.coord_equal()
        + pn.coord_flip()
    )
    display(plot)
    
    path = snakemake.params["output_basedir"] + f"/r2_bar_plot.{feature_x}__vs__{feature_y}"
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)


# %% [markdown]
# ## barplot num. significants

# %% [raw]
# plot_df = num_significant_associations
# plot_df = plot_df.assign(
#     covariates=plot_df["covariates"].str.replace("_", "\n+ "),
#     phenotype_col=(
#         plot_df["phenotype_col"]
#         .str.replace(r"_(f\d+_.*)", r"\n(\1)", regex=True)
#         .str.split("\n")
#         # .str.replace(r"_", r" ", regex=True)
#         .apply(
#             lambda s: "\n".join([
#                 textwrap.fill(s[0].replace("_", " "), 12, break_long_words=False),
#                 *s[1:]
#             ])
#         )
#         .astype("string[pyarrow]")
#     ),
# )
#
# # crop_pvalue = 10 ** -10

# %%
plot_df = quantiles_pd_df
plot_df = plot_df.astype({
    "prediction_quantile": "str",
    "measurement_quantile": "str",
})
plot_df = plot_df.assign(**{
    "difference_to_restricted_model": plot_df["full_model"] - plot_df["restricted_model"],
    "proportional_difference_to_restricted_model": (plot_df["full_model"] / plot_df["total"]) - (plot_df["restricted_model"] / plot_df["total"]),
})

# %%
plot_df.columns

# %%
grouping = ['measurement_quantile', 'prediction_quantile', 'phenotype_col', 'feature_set', 'covariates', 'bound', 'total']
keys = plot_df["feature_set"].unique().tolist()

unstacked_plot_df = plot_df.set_index(grouping)['full_model'].unstack("feature_set").reset_index(level="total")
unstacked_plot_df

# %%
import itertools

# list(itertools.combinations(keys, 2))
for feature_x, feature_y in list(itertools.product(keys, keys)):
    if feature_x == feature_y:
        continue
    import mizani
    
    subset_plot_df = unstacked_plot_df.assign(**{
        f"difference_to_{feature_x}": unstacked_plot_df[feature_y] - unstacked_plot_df[feature_x],
        f"proportional_difference_to_{feature_x}": (unstacked_plot_df[feature_y] / unstacked_plot_df["total"]) - (unstacked_plot_df[feature_x] / unstacked_plot_df["total"]),
    })
    subset_plot_df = subset_plot_df.reset_index()
    subset_plot_df = subset_plot_df.assign(
        covariates=subset_plot_df["covariates"].str.replace("_", "\n+ "),
        phenotype_col=(
            subset_plot_df["phenotype_col"]
            .str.replace(r"_(f\d+_.*)", r"\n(\1)", regex=True)
            .str.split("\n")
            # .str.replace(r"_", r" ", regex=True)
            .apply(
                lambda s: "\n".join([
                    textwrap.fill(s[0].replace("_", " "), 12, break_long_words=False),
                    *s[1:]
                ])
            )
            .astype("string[pyarrow]")
        ),
    )

    plot = (
        pn.ggplot(subset_plot_df.reset_index(), pn.aes(x="prediction_quantile", y="measurement_quantile", fill=f"proportional_difference_to_{feature_x}"))
        + pn.geom_tile(
            # pn.aes(width=.95, height=.95)
        )
        + pn.scale_fill_gradient2(
            # low = muted("red"),
            # mid = "white",
            # high = muted("blue"),
            midpoint = 0,
            labels=mizani.formatters.percent_format()
        )
        + pn.labs(
            fill=f"proportional difference",
            title=f"difference of individuals at risk for '{feature_y}'\n(baseline: Age+Sex+PRS+{feature_x})"
        )
        + pn.theme(
            legend_text=pn.element_text(linespacing=1.4),
            figure_size=(8, 100),
            axis_text_x=pn.element_text(
                rotation=45,
                hjust=1
            ),
            # strip_text_y=pn.element_text(
            #     rotation=0,
            # ),
            title=pn.element_text(linespacing=1.4, vjust=-10),
        )
        + pn.facet_grid(
            "phenotype_col ~ bound",
            # scales="free_y"
            # scales="free"
        )
        + pn.coord_equal()
        + pn.coord_flip()
    )
    display(plot)
    
    path = snakemake.params["output_basedir"] + f"/diff_individuals_at_risk.heatmap.{feature_x}__vs__{feature_y}"
    pn.ggsave(plot, path + ".png", dpi=DPI, limitsize=False)
    pn.ggsave(plot, path + ".pdf", dpi=DPI, limitsize=False)


# %%
import itertools

# list(itertools.combinations(keys, 2))
for feature_x, feature_y in list(itertools.product(keys, keys)):
    if feature_x == feature_y:
        continue
    
    subset_plot_df = unstacked_plot_df.query("prediction_quantile == measurement_quantile").reset_index()
    subset_plot_df = subset_plot_df.assign(**{
        f"difference_to_{feature_x}": subset_plot_df[feature_y] - subset_plot_df[feature_x],
        f"proportional_difference_to_{feature_x}": (subset_plot_df[feature_y] / subset_plot_df["total"]) - (subset_plot_df[feature_x] / subset_plot_df["total"]),
        f"true_positive_rate_{feature_x}": (subset_plot_df[feature_x] / subset_plot_df["total"]),
        f"true_positive_rate_{feature_y}": (subset_plot_df[feature_y] / subset_plot_df["total"]),
        "quantile": subset_plot_df["measurement_quantile"].astype("str") + " (" + subset_plot_df["bound"] + ")",
        "phenotype_col": subset_plot_df["phenotype_col"].str.replace("_", " "),
    })
    subset_plot_df = subset_plot_df.reset_index()
    subset_plot_df = subset_plot_df.assign(
        covariates=subset_plot_df["covariates"].str.replace("_", "\n+ "),
        phenotype_col=(
            subset_plot_df["phenotype_col"]
            .str.replace(r"_(f\d+_.*)", r"\n(\1)", regex=True)
            .str.split("\n")
            # .str.replace(r"_", r" ", regex=True)
            .apply(
                lambda s: "\n".join([
                    textwrap.fill(s[0].replace("_", " "), 12, break_long_words=False),
                    *s[1:]
                ])
            )
            .astype("string[pyarrow]")
        ),
    )

    plot = (
        pn.ggplot(subset_plot_df, pn.aes(
            x=f"true_positive_rate_{feature_x}",
            y=f"true_positive_rate_{feature_y}",
            fill=f"phenotype_col",
            shape="quantile",
        ))
        + pn.geom_point(
            # pn.aes(width=.95, height=.95),
            size=3,
        )
        + pn.geom_abline(slope=1, linetype="dashed")
        # + pn.scale_fill_gradient2(
        #     # low = muted("red"),
        #     # mid = "white",
        #     # high = muted("blue"),
        #     midpoint = 0,
        #     labels=mizani.formatters.percent_format()
        # )
        + pn.labs(
            x=f"""true positive rate ({feature_x.replace("_", " ")})""",
            y=f"""true positive rate ({feature_y.replace("_", " ")})""",
            # title=f"difference of individuals at risk for '{feature_y}'\n(baseline: Age+Sex+PRS+{feature_x})"
        )
        + pn.theme(
            legend_text=pn.element_text(linespacing=1.4),
            # figure_size=(8, 8),
            # axis_text_x=pn.element_text(
            #     rotation=45,
            #     hjust=1
            # ),
            # strip_text_y=pn.element_text(
            #     rotation=0,
            # ),
            title=pn.element_text(linespacing=1.4, vjust=-10),
        )
        + pn.coord_equal()
        # + pn.coord_flip()
    )
    display(plot)
    
    path = snakemake.params["output_basedir"] + f"/true_positive_rate.scatter.{feature_x}__vs__{feature_y}"
    pn.ggsave(plot, path + ".png", dpi=DPI, limitsize=False)
    pn.ggsave(plot, path + ".pdf", dpi=DPI, limitsize=False)


# %%
snakemake.output

# %%
