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

import os
import numpy as np
import pandas as pd

import json
import yaml

import pyspark
import pyspark.sql.types as t
import pyspark.sql.functions as f

import polars as pl

import re
import patsy

import plotnine as pn

import matplotlib
import matplotlib.pyplot as plt

import lightgbm as lgb

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

from rep.notebook_init import setup_plot_style
setup_plot_style()

# %%
# we need ray to efficiently share the covariates df in memory
# use_ray = True
use_ray = False

# %%
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
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__polygenic_risk_score',
        default_wildcards={
            "phenotype_col": "c_reactive_protein",
            # "phenotype_col": "standing_height",
            # "phenotype_col": "glycated_haemoglobin_hba1c",
            # "phenotype_col": "Lipoprotein_A",
            # "phenotype_col": "BodyMassIndex",
            # "phenotype_col": "Triglycerides",
            # "phenotype_col": "LDL_direct",
            # "phenotype_col": "systolic_blood_pressure",
            # "phenotype_col": "HDL_cholesterol",
            # "feature_set": "LOFTEE_pLoF",
            # "feature_set": "AbExp_all_tissues",
            "feature_set": "LOFTEE_pLoF",
            "covariates": "sex_age_genPC_CLMP_PRS",
            # "covariates": "sex_age_genPC_CLMP",
            # "covariates": "sex_age_genPC",
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

# %% [markdown]
# # Load configuration

# %%
with open(snakemake.input["featureset_config"], "r") as fd:
    config = yaml.safe_load(fd)

# %%
snakemake.input["featureset_config"]

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

# %%
covariates_df_columns = covariates_df.columns
#covariates_df_columns

# %% [raw]
# ## clumping

# %% [raw]
# clumping_variants_df = (
#     pl.read_parquet(snakemake.input["clumping_variants_pq"])
#     .filter(pl.col("gene").is_in(significant_genes))
# )
# clumping_variants_df

# %% [raw]
# clumping_variants = clumping_variants_df["variant"].unique().to_list()

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

# %% [markdown]
# ### pivot genes

# %%
import itertools

# %%
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
#pivot_features_df.printSchema()

# %%
gene_features = [f"{g}_{s}" for g, s in itertools.product(significant_genes, features)]

# %% [markdown]
# ### join everything

# %%
full_features_df = (
    covariates_df
    .join(pivot_features_df, on="individual", how="left")
    .fillna(0, subset=pivot_features_df.columns)
)
#full_features_df.printSchema()

# %% [markdown]
# # select variables

# %%
def format_formula(formula, add_clumping=False, clumping_variants=None):
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
formatted_restricted_formula = patsy.ModelDesc.from_formula(restricted_formula)
formatted_restricted_formula

# %%
formatted_full_formula = format_formula(
    formula=restricted_formula,
    add_clumping=False, # do not add clumping variants!
)
if len(gene_features) > 0:
    formatted_full_formula.rhs_termlist += [
        patsy.Term([patsy.EvalFactor(f"Q('{c}')")]) for c in gene_features
    ]

# %%
formatted_full_formula.rhs_termlist[-10:]

# %%
import patsy

def get_variables_from_formula(formula, lhs=True, rhs=True):
    if not isinstance(formula, patsy.ModelDesc):
        model_desc = patsy.ModelDesc.from_formula(formula)
    else:
        model_desc = formula

    covariate_cols = [
        *([factor.name() for term in model_desc.lhs_termlist for factor in term.factors] if lhs else []),
        *([factor.name() for term in model_desc.rhs_termlist for factor in term.factors] if rhs else []),
    ]
    
    # find Q($i) -> $i
    q_regex = re.compile(r"""Q\(['"](.*)['"]\)""")
    
    parsed_covariate_cols = []
    for cov in covariate_cols:
        q_matches = q_regex.findall(cov)
        
        if len(q_matches) > 0:
            parsed_covariate_cols += q_matches
        else:
            parsed_covariate_cols.append(cov)

    # deduplicate
    parsed_covariate_cols = list(dict.fromkeys(parsed_covariate_cols))
    
    return parsed_covariate_cols


# %%
full_variables = get_variables_from_formula(formatted_full_formula, lhs=False)
full_variables

# %%
restricted_variables = get_variables_from_formula(formatted_restricted_formula, lhs=False)
restricted_variables

# %%
# Model with AGE+SEX+PC
basic_variables = [c for c in restricted_variables if not c.startswith("PGS")]
basic_variables

# %% [markdown]
# # write features

# %%
prs_features_df = (
    full_features_df
    .select([
        "individual",
        phenotype_col,
        *full_variables,
    ])
    .sort("individual")
)

# %%
snakemake.output

# %%
#prs_features_df.write.parquet(snakemake.output["prs_features_pq"], mode="overwrite")

# %%
#prs_features_df = spark.read.parquet(snakemake.output["prs_features_pq"])

# %%
prs_features_pd_df = (
    prs_features_df
    .withColumn("sex_f31_0_0", f.col("sex_f31_0_0").cast(t.ByteType()))
    .toPandas()
)
prs_features_pd_df

# %%
spark.sparkContext._jvm.System.gc()

# %% [markdown]
# # read samples

# %%
train_test_split = pd.read_parquet(snakemake.input["train_test_split_pq"])
train_test_split

# %%
samples = train_test_split[["individual", "fold"]]
samples

# %% [markdown]
# ## join with prs features

# %%
prs_features_pd_df = prs_features_pd_df.merge(samples, on="individual", how="inner")
prs_features_pd_df

# %%
prs_features_pd_df = prs_features_pd_df.dropna()
prs_features_pd_df

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
train_samples = samples.loc[samples["fold"] != "association_testing", "individual"].values
train_samples

# %% [raw]
# test_samples = samples.loc[samples["fold"] == "association_testing", "individual"].values
# test_samples

# %%
train_df = prs_features_pd_df[prs_features_pd_df["individual"].isin(train_samples)]
train_df

# %%
cv_split = sklearn.model_selection.PredefinedSplit(train_df["fold"].str.split(" ").apply(lambda s: s[-1]).astype("int"))
cv_split

# %%
cv_split.get_n_splits()

# %%
full_variables

# %%
cv_pred_dfs = []
for fold, (train_fold_index, test_fold_index) in enumerate(cv_split.split()):
    print(f"training fold {fold}...")
    X_fold_train = train_df.iloc[train_fold_index]
    y_fold_train = train_df.iloc[train_fold_index][phenotype_col]
    X_fold_test = train_df.iloc[test_fold_index][full_variables]
    # y_fold_test = train_df.iloc[test_fold_index][phenotype_col]
    
    gbm_full = lgb.LGBMRegressor()
    gbm_full.fit(X_fold_train[full_variables], y_fold_train)
    full_pred_test = gbm_full.predict(X_fold_test)
    
    gbm_restricted = lgb.LGBMRegressor()
    gbm_restricted.fit(X_fold_train[restricted_variables], y_fold_train)
    restricted_pred_test = gbm_restricted.predict(X_fold_test[restricted_variables])
    
    gbm_basic = lgb.LGBMRegressor()
    gbm_basic.fit(X_fold_train[basic_variables], y_fold_train)
    basic_pred_test = gbm_basic.predict(X_fold_test[basic_variables])
    
    res = train_df.iloc[test_fold_index][["individual", "fold", phenotype_col]].assign(**{
        "full_model_pred": full_pred_test,
        "restricted_model_pred": restricted_pred_test,
        "basic_model_pred": basic_pred_test,
    })
    cv_pred_dfs.append(res)

# %%
cv_pred_df = (
    pd.concat(cv_pred_dfs)
    .rename(columns={phenotype_col: "measurement"})
    .assign(**snakemake.wildcards)
)
cv_pred_df

# %%
# Save predictions
(
    cv_pred_df
    .to_parquet(snakemake.output["predictions_pq"], index=False)
)


# %%
def get_r2_scores(df, models=["full", "restricted", "basic"]):
    r2_scores = {}
    for model in models:
        y_test = df["measurement"]
        pred = df[f"{model}_model_pred"]

        r2 = sklearn.metrics.r2_score(y_test, pred)
        r2_scores[model] = r2
    return r2_scores


# %%
get_r2_scores(cv_pred_df)

# %%
fold_r2_scores = (
    cv_pred_df
    .groupby("fold")
    .apply(lambda df: pd.Series(get_r2_scores(df)))
    .rename(columns={
        "full": "full_model_r2",
        "restricted": "restricted_model_r2",
        "basic": "basic_model_r2",
    })
    .assign(**snakemake.wildcards)
    .reset_index()
)
fold_r2_scores

# %%
# Save r2 scores
fold_r2_scores.to_parquet(snakemake.output["r2_scores_pq"], index=False)
fold_r2_scores.to_csv(snakemake.output["r2_scores_tsv"], sep="\t", index=False)

# %%
# Calc PRC
prc_df = []
for model in ["basic", "restricted", "full"]:
    for extreme in ["bottom", "top"]:
        for percentile in [0.01 , 0.05, 0.1, 0.2]:
            if extreme == "top":
                measurement_quantile = cv_pred_df["measurement"].quantile(1-percentile)
                prc = sklearn.metrics.precision_recall_curve(cv_pred_df["measurement"].ge(measurement_quantile), cv_pred_df[f"{model}_model_pred"])
            else:
                measurement_quantile = cv_pred_df["measurement"].quantile(percentile)
                prc = sklearn.metrics.precision_recall_curve(~cv_pred_df["measurement"].ge(measurement_quantile), (-1) * cv_pred_df[f"{model}_model_pred"])
            prc_df.append(pd.DataFrame({"extreme": extreme, "percentile": percentile, "method": model, "precision" : prc[0], "recall" : prc[1], "auPRC": sklearn.metrics.auc(prc[1], prc[0])}))
prc_df = pd.concat(prc_df)
prc_df["method"] = prc_df["method"].replace({"basic": "Age+Sex+PC", "restricted": "Age+Sex+PC+PRS", "full": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}"})
prc_df = prc_df.assign(**snakemake.wildcards)

# %%
# Save prc 
prc_df.query("method=='Age+Sex+PC' or method=='Age+Sex+PC+PRS'").to_parquet(snakemake.output["precision_recall_baseline_pq"], index=False)
prc_df.query(f"method=='Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}'").to_parquet(snakemake.output["precision_recall_full_pq"], index=False)

# %%
prc_df

# %% [raw]
#
