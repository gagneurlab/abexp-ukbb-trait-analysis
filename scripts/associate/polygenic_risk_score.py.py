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

# %% {"tags": []}
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

# %% {"tags": []}
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

from rep.notebook_init import setup_plot_style
setup_plot_style()

# %% {"tags": []}
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

# %% {"tags": []}
ray_context

# %% {"tags": []}
spark

# %% {"tags": []}
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %% {"tags": []}
#del snakemake

# %% {"tags": []}
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

# %% {"tags": []}
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %% {"tags": []}
if "plot_dpi" in snakemake.params:
    DPI = snakemake.params["plot_dpi"]
else:
    DPI=450

# %% [markdown] {"tags": []}
# # Load configuration

# %% {"tags": []}
with open(snakemake.input["featureset_config"], "r") as fd:
    config = yaml.safe_load(fd)

# %% {"tags": []}
snakemake.input["featureset_config"]

# %% {"tags": []}
print(json.dumps(config, indent=2, default=str))

# %% {"tags": []}
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %% {"tags": []}
phenocode = config["covariates"]["phenocode"]
phenocode

# %% {"tags": []}
restricted_formula = config["covariates"]["restricted_formula"]
print(restricted_formula)

# %% {"tags": []}
pvalue_cutoff = 0.05

# %% [markdown]
# # Read association results

# %% {"tags": []}
snakemake.input["associations_pq"]

# %% {"tags": []}
regression_results_df = spark.read.parquet(snakemake.input["associations_pq"])

# %% {"tags": []}
correction_value = regression_results_df.count()
print(f"Correcting for {correction_value} association tests...")

# %% {"tags": []}
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

# %% {"tags": []}
regression_results_df.printSchema()

# %% [markdown]
# ## extract significant genes

# %% {"tags": []}
significant_genes = (
    regression_results_df
    .filter(f.col("padj") < f.lit(pvalue_cutoff))
    .select("gene")
    .toPandas()["gene"].tolist()
)
significant_genes

# %% [markdown]
# # Read covariates

# %% {"tags": []}
covariates_df = spark.read.parquet(snakemake.input["covariates_pq"])

# %% {"tags": []}
covariates_df_columns = covariates_df.columns
#covariates_df_columns

# %% [raw] {"tags": []}
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

# %% {"tags": []}
feature_dfs = {}
for feature_name, path in config["feature_sets"].items():
    feature_dfs[feature_name] = (
        spark.read.parquet(path + "/data.parquet")
        .filter(f.col("individual").isNotNull())
        .filter(f.col("gene").isNotNull())
    )

# %% {"tags": []}
len(feature_dfs)

# %% {"tags": []}
from rep.data import join_featuresets

# %% {"tags": []}
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

# %% {"tags": []}
# features_df.select("gene").distinct().count()

# %% {"tags": []}
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

# %% {"tags": []}
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
#pivot_features_df.printSchema()

# %% {"tags": []}
gene_features = [f"{g}_{s}" for g, s in itertools.product(significant_genes, features)]

# %% [markdown]
# ### join everything

# %% {"tags": []}
full_features_df = (
    covariates_df
    .join(pivot_features_df, on="individual", how="left")
    .fillna(0, subset=pivot_features_df.columns)
)
#full_features_df.printSchema()

# %% [markdown]
# # select variables

# %% {"tags": []}
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


# %% {"tags": []}
formatted_restricted_formula = patsy.ModelDesc.from_formula(restricted_formula)
formatted_restricted_formula

# %% {"tags": []}
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


# %% {"tags": []}
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
full_variables = get_variables_from_formula(formatted_full_formula, lhs=False)
full_variables

# %% {"tags": []}
restricted_variables = get_variables_from_formula(formatted_restricted_formula, lhs=False)
restricted_variables

# %% {"tags": []}
# Model with AGE+SEX+PC
basic_variables = [c for c in restricted_variables if not c.startswith("PGS")]
basic_variables

# %% [markdown]
# # write features

# %% {"tags": []}
prs_features_df = (
    full_features_df
    .select([
        "individual",
        phenotype_col,
        *full_variables,
    ])
    .sort("individual")
)

# %% {"tags": []}
snakemake.output

# %% {"tags": []}
#prs_features_df.write.parquet(snakemake.output["prs_features_pq"], mode="overwrite")

# %% {"tags": []}
#prs_features_df = spark.read.parquet(snakemake.output["prs_features_pq"])

# %% {"tags": []}
prs_features_pd_df = (
    prs_features_df
    .withColumn("sex_f31_0_0", f.col("sex_f31_0_0").cast(t.ByteType()))
    .toPandas()
)
prs_features_pd_df

# %% {"tags": []}
spark.sparkContext._jvm.System.gc()

# %% [markdown]
# # read samples

# %% {"tags": []}
train_test_split = pd.read_parquet(snakemake.input["train_test_split_pq"])
train_test_split

# %% {"tags": []}
samples = train_test_split[["individual", "fold"]]
samples

# %% [markdown]
# ## join with prs features

# %% {"tags": []}
prs_features_pd_df = prs_features_pd_df.merge(samples, on="individual", how="inner")
prs_features_pd_df

# %% {"tags": []}
prs_features_pd_df = prs_features_pd_df.dropna()
prs_features_pd_df

# %% [markdown]
# # perform regression

# %% {"tags": []}
import sklearn
import sklearn.pipeline
import sklearn.linear_model
import sklearn.preprocessing
import sklearn.metrics

# %% [markdown]
# ## prepare training data

# %% {"tags": []}
train_samples = samples.loc[samples["fold"] != "association_testing", "individual"].values
train_samples

# %% [raw] {"tags": []}
# test_samples = samples.loc[samples["fold"] == "association_testing", "individual"].values
# test_samples

# %% {"tags": []}
train_df = prs_features_pd_df[prs_features_pd_df["individual"].isin(train_samples)]
train_df

# %% {"tags": []}
cv_split = sklearn.model_selection.PredefinedSplit(train_df["fold"].str.split(" ").apply(lambda s: s[-1]).astype("int"))
cv_split

# %% {"tags": []}
cv_split.get_n_splits()

# %% {"tags": []}
full_variables

# %% {"tags": []}
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

# %% {"tags": []}
cv_pred_df = (
    pd.concat(cv_pred_dfs)
    .rename(columns={phenotype_col: "measurement"})
    .assign(**snakemake.wildcards)
)
cv_pred_df

# %% {"tags": []}
# Save predictions
(
    cv_pred_df
    .to_parquet(snakemake.output["predictions_pq"], index=False)
)


# %% {"tags": []}
def get_r2_scores(df, models=["full", "restricted", "basic"]):
    r2_scores = {}
    for model in models:
        y_test = df["measurement"]
        pred = df[f"{model}_model_pred"]

        r2 = sklearn.metrics.r2_score(y_test, pred)
        r2_scores[model] = r2
    return r2_scores


# %% {"tags": []}
get_r2_scores(cv_pred_df)

# %% {"tags": []}
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

# %% {"tags": []}
# Save r2 scores
fold_r2_scores.to_parquet(snakemake.output["r2_scores_pq"], index=False)
fold_r2_scores.to_csv(snakemake.output["r2_scores_tsv"], sep="\t", index=False)

# %% {"tags": []}
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

# %% {"tags": []}
# Save prc 
prc_df.query("method=='Age+Sex+PC' or method=='Age+Sex+PC+PRS'").to_parquet(snakemake.output["precision_recall_baseline_pq"], index=False)
prc_df.query(f"method=='Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}'").to_parquet(snakemake.output["precision_recall_full_pq"], index=False)

# %% {"tags": []}
prc_df

# %% [raw]
#
