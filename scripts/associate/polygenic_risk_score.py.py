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
import polars as pl

# %%
import re
import patsy

# %%
import plotnine as pn
# import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

# %%
from rep.notebook_init import setup_plot_style
setup_plot_style()

# %%
# import os
# # os.environ["RAY_ADDRESS"] = os.environ.get("RAY_ADDRESS", 'ray://192.168.16.30:10001')
# os.environ["RAY_ADDRESS"] = 'ray://192.168.16.28:10001'
# os.environ["RAY_ADDRESS"]

# %% [raw] {"tags": []}
#

# %%
# we need ray to efficiently share the covariates df in memory
# use_ray = True
use_ray = False

# %% {"tags": []}
ray = None
ray_context = None
if use_ray:
    import ray
    from rep.notebook_init import init_ray, init_spark_on_ray
    
    ray_context = init_ray(
        plasma_store_memory_fraction=0.3
    )

    spark = init_spark_on_ray(
        # executor_cores=128,
        executor_memory_overhead=0.9,
        # configs={
        #     "spark.default.parallelism": int(ray.cluster_resources()["CPU"]),
        #     "spark.sql.shuffle.partitions": 2048,
        # }
        enable_glow=False,
    )
else:
    from rep.notebook_init import init_spark
    spark = init_spark(enable_glow=False)

# %%
ray_context

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
        rule_name = 'associate__polygenic_risk_score',
        default_wildcards={
            # "phenotype_col": "triglycerides_f30870_0_0",
            # "phenotype_col": "standing_height_f50_0_0",
            # "phenotype_col": "body_mass_index_bmi_f21001_0_0",
            # "phenotype_col": "systolic_blood_pressure_automated_reading_f4080_0_0",
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            # "feature_set": "LOFTEE_pLoF",
            "feature_set": "AbExp_all_tissues",
            # "feature_set": "LOFTEE_pLoF",
            "covariates": "sex_age_genPC_CLMP_PRS",
            # "covariates": "sex_age_genPC_CLMP",
            # "covariates": "sex_age_genPC",
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

# %%
pvalue_cutoff = 0.05

# %% [markdown]
# # Read association results

# %%
snakemake.input["associations_pq"]

# %%
regression_results_df = spark.read.parquet(snakemake.input["associations_pq"])

# %%
correction_value = regression_results_df.count()
print(f"Correcting for {correction_value} association tests...")

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

# %% [markdown]
# ## extract significant genes

# %%
significant_genes = (
    regression_results_df
    .filter(f.col("padj") < f.lit(pvalue_cutoff))
    .select("gene")
    .toPandas()["gene"].tolist()
)
significant_genes

# %% [markdown]
# # Read covariates

# %%
covariates_df = spark.read.parquet(snakemake.input["covariates_pq"])

# %% {"tags": []}
covariates_df_columns = covariates_df.columns
covariates_df_columns

# %% [markdown] {"tags": []}
# ## clumping

# %%
clumping_variants_df = (
    pl.read_parquet(snakemake.input["clumping_variants_pq"])
    .filter(pl.col("gene").is_in(significant_genes))
)
clumping_variants_df

# %%
clumping_variants = clumping_variants_df["variant"].unique().to_list()

# %% [markdown]
# ## feature dfs

# %%
feature_dfs = {}
for feature_name, path in config["feature_sets"].items():
    feature_dfs[feature_name] = (
        spark.read.parquet(path + "/data.parquet")
        .filter(f.col("individual").isNotNull())
        .filter(f.col("gene").isNotNull())
    )

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
    index_cols=["individual", "gene"],
    fill_values=fill_values,
    join="outer",
)
features_df = (
    features_df
    .filter(f.col("gene").isin(significant_genes))
)
features_df.printSchema()

# %%
# features_df.select("gene").distinct().count()

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
# ### pivot genes

# %%
import itertools

# %% {"tags": []}
pivot_features_df = (
    renamed_features_df
    .groupby("individual")
    .pivot("gene", values=significant_genes)
    .agg(
        # f.first(f.col("y_pred_proba")).alias("AbExp")
        f.struct(
            *[
                f.coalesce(
                    f.first(f.col(c)),
                    f.lit(0)
                ).alias(c) for c in features
            ],
        )
    )
    .select(
        "individual",
        # rename all columns to url-encoded versions in order to escape special characters
        *[f.col(g)[s].alias(f"{g}_{s}") for g, s in itertools.product(significant_genes, features)]
    )
)
pivot_features_df.printSchema()

# %%
gene_features = [f"{g}_{s}" for g, s in itertools.product(significant_genes, features)]

# %% [markdown]
# ### join everything

# %% {"tags": []}
full_features_df = (
    covariates_df
    .join(pivot_features_df, on="individual", how="left")
    .fillna(0, subset=pivot_features_df.columns)
)
full_features_df.printSchema()


# %% [markdown]
# # select variables

# %%
def format_formula(formula, add_clumping=True, clumping_variants=clumping_variants):
    if not isinstance(formula, patsy.ModelDesc):
        model_desc = patsy.ModelDesc.from_formula(formula)
    else:
        model_desc = formula
    
    if add_clumping:
        if len(clumping_variants) > 0:
            model_desc.rhs_termlist += [
                patsy.Term([patsy.EvalFactor(f"Q('{c}')")]) for c in clumping_variants
            ]
    
    return model_desc


# %%
# formatted_restricted_formula = patsy.ModelDesc.from_formula(restricted_formula)

# %%
formatted_full_formula = format_formula(
    formula=restricted_formula,
    add_clumping=False, # do not add clumping variants!
)
if len(gene_features) > 0:
    formatted_full_formula.rhs_termlist += [
        patsy.Term([patsy.EvalFactor(f"Q('{c}')")]) for c in gene_features
    ]

# %% {"tags": []}
formatted_full_formula.rhs_termlist[-10:]


# %%
def get_variables_from_formula(formula, lhs=True, rhs=True):
    if not isinstance(formula, patsy.ModelDesc):
        model_desc = patsy.ModelDesc.from_formula(formula)
    else:
        model_desc = formula

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


# %% {"tags": []}
full_variables = get_variables_from_formula(formatted_full_formula)
len(full_variables)

# %% [markdown]
# # write features

# %%
prs_features_df = (
    full_features_df
    .select([
        "individual",
        *full_variables,
    ])
    .sort("individual")
)

# %%
snakemake.output

# %%
prs_features_df.write.parquet(snakemake.output["prs_features_pq"], mode="overwrite")

# %%
prs_features_df = spark.read.parquet(snakemake.output["prs_features_pq"])

# %%
prs_features_pd_df = prs_features_df.toPandas()

# %%
spark.sparkContext._jvm.System.gc()

# %%
prs_features_pd_df.head()

# %%

# %% [markdown]
# # read samples

# %%
samples = pd.read_parquet(snakemake.input["samples_pq"])
samples

# %%
train_samples = samples.loc[samples["fold"] != "test", "individual"].values
train_samples

# %%
test_samples = samples.loc[samples["fold"] == "test", "individual"].values
test_samples

# %% [markdown]
# ## join with prs features

# %%
prs_features_pd_df = prs_features_pd_df.merge(samples, on="individual", how="inner")
prs_features_pd_df

# %% [markdown]
# ## prepare splits

# %%
non_test_idx

# %% [markdown]
# # perform regression

# %%
import sklearn
import sklearn.pipeline
import sklearn.linear_model
import sklearn.preprocessing
import sklearn.metrics

# %% [markdown]
# ## prepare training data

# %%
train_df = prs_features_pd_df[prs_features_pd_df["individual"].isin(train_samples)]
train_df

# %%
cv_split = sklearn.model_selection.PredefinedSplit(train_df["fold"].str.split(" ").apply(lambda s: s[-1]).astype("int"))
cv_split

# %%
cv_split.get_n_splits()

# %%
test_df = prs_features_pd_df[prs_features_pd_df["individual"].isin(test_samples)]
test_df

# %%
X_train = train_df[restricted_variables].drop(columns=phenotype_col)
y_train = train_df[phenotype_col]

# %%
X_test = test_df[restricted_variables].drop(columns=phenotype_col)
y_test = test_df[phenotype_col]

# %% [markdown]
# ## train models

# %%
# TODO for Jonas

# %% [raw]
# ### restricted model

# %% [raw]
# model_pipeline = sklearn.pipeline.Pipeline([
#     ("normalize", sklearn.preprocessing.StandardScaler()),
#     ("linear_model", sklearn.linear_model.ElasticNetCV(cv=cv_split)),
# ])

# %% [raw]
# restricted_model = model_pipeline.fit(X_train, y_train)
# restricted_model

# %% [raw]
# restricted_pred_test = restricted_model.predict(X_test)
# restricted_pred_test

# %% [raw]
# sklearn.metrics.r2_score(y_test, restricted_pred_test)

# %% [raw]
# ### full model

# %%
model_pipeline = sklearn.pipeline.Pipeline([
    ("normalize", sklearn.preprocessing.StandardScaler()),
    ("linear_model", sklearn.linear_model.ElasticNetCV(cv=cv_split)),
])

# %%
full_model = model_pipeline.fit(X_train, y_train)
full_model

# %%
full_pred_test = full_model.predict(X_test)
full_pred_test

# %%
sklearn.metrics.r2_score(y_test, full_pred_test)

# %%

# %%

# %% {"tags": []}
import statsmodels.formula.api as smf

# %%
restricted_model = smf.ols(
    formatted_restricted_formula,
    data = prs_features_pd_df
).fit()

# %%
full_model = smf.ols(
    formatted_full_formula,
    data = prs_features_pd_df
).fit()

# %%
lr_stat, lr_pval, lr_df_diff = full_model.compare_lr_test(restricted_model)

# %%
lr_pval

# %%
lr_df_diff

# %%
prs_features_pd_df.dropna().shape

# %%
prs_features_pd_df.shape

# %%
pred_df = (
    prs_features_pd_df[["individual", phenotype_col]]
    .rename(columns={phenotype_col: "measurement"})
    .assign(**{
        "phenotype_col": phenotype_col,
        "full_model_pred": full_model.predict(prs_features_pd_df),
        "restricted_model_pred": restricted_model.predict(prs_features_pd_df),
    })
)
pred_df

# %%
pred_df = (
    pred_df
    .assign(**{
        "residuals": pred_df["measurement"] - pred_df["restricted_model_pred"],
        "full - restricted": pred_df["full_model_pred"] - pred_df["restricted_model_pred"],
        
    })
)
pred_df

# %%
restricted_model.rsquared

# %%
full_model.rsquared

# %% [markdown]
# # save model + predictions

# %%
snakemake.output

# %% [markdown]
# ## predictions

# %%
pred_df.to_parquet(snakemake.output["predictions_pq"], index=False)

# %% [markdown]
# ## summary

# %%
restricted_model_summary = restricted_model.summary2()
full_model_summary = full_model.summary2()

# %%
with open(snakemake.output["restricted_summary_txt"], "w") as fd:
    fd.write(restricted_model_summary.as_text())
    fd.write("\n")

# %%
with open(snakemake.output["full_summary_txt"], "w") as fd:
    fd.write(full_model_summary.as_text())
    fd.write("\n")

# %% [markdown]
# ## params

# %%
restricted_params_df = (
    restricted_model_summary.tables[1]
    .rename_axis(index="term")
    .reset_index()
)
restricted_params_df

# %%
restricted_params_df.to_parquet(snakemake.output["restricted_params_pq"], index=False)

# %%
full_params_df = (
    full_model_summary.tables[1]
    .rename_axis(index="term")
    .reset_index()
)
full_params_df

# %%
full_params_df.to_parquet(snakemake.output["full_params_pq"], index=False)

# %% [markdown]
# # plotting

# %%
# TODO for Jonas

# %% [markdown]
# ## phenotype correlation

# %%
plot_df = pred_df[["phenotype_col", "measurement", "full_model_pred", "restricted_model_pred"]]
plot_df = plot_df.melt(id_vars=["phenotype_col", "measurement"])

plot = (
    pn.ggplot(plot_df, pn.aes(y="measurement", x="value"))
    + pn.geom_bin_2d(bins=100)
    + pn.geom_smooth(method="lm", color="red")
    + pn.facet_grid("phenotype_col ~ variable")
    + pn.scale_fill_continuous(trans = "log10")
    + pn.coord_equal()
    + pn.labs(
        x="prediction",
    )
)
display(plot)

# %%
snakemake.params["output_basedir"]

# %%
path = snakemake.params["output_basedir"] + "/phenotype_correlation"
pn.ggsave(plot, path + ".png", dpi=DPI)
pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [raw]
# from IPython.display import Image
# display(Image(snakemake.params["output_basedir"] + "/phenotype_correlation.png"))

# %% [markdown]
# ## residual correlation

# %%
snakemake.wildcards

# %%
import sklearn.metrics as skmetrics
import textwrap

plot_r2 = skmetrics.r2_score(pred_df["residuals"], pred_df["full - restricted"])

plot = (
    pn.ggplot(pred_df, pn.aes(y="residuals", x="full - restricted"))
    + pn.geom_bin_2d(bins=100)
    + pn.geom_smooth(method="lm", color="red")
    + pn.scale_fill_continuous(trans = "log10")
    + pn.labs(
        x="prediction\n(full model - restricted model)",
        y=f"residuals\n(phenotype - restricted model)",
        title=textwrap.dedent(
            f"""
            Polygenic score predicting '{phenotype_col}'
            Covariates: '{snakemake.wildcards['covariates'].replace('_', '+')}'
            r2: '{"%.4f" % plot_r2}'
            (Every point represents one individual)\
            """
        ).strip(),
        fill='Nr. of individuals'
    )
    + pn.theme(title=pn.element_text(linespacing=1.4))
    + pn.coord_equal()
)

# %%
if len(restricted_variables) != len(full_variables):
    display(plot)

# %%
snakemake.params["output_basedir"]

# %%
path = snakemake.params["output_basedir"] + "/residual_correlation"
if len(restricted_variables) != len(full_variables):
    pn.ggsave(plot, path + ".png", dpi=DPI)
    pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %%
