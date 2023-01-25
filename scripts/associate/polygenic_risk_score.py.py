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
            "phenotype_col": "standing_height",
            #"phenotype_col": "glycated_haemoglobin_hba1c",
            #"phenotype_col": "Lipoprotein_A",
            #"phenotype_col": "BodyMassIndex",
            #"phenotype_col": "Triglycerides",
            #"phenotype_col": "LDL_direct",
            #"phenotype_col": "systolic_blood_pressure",
            #"phenotype_col": "HDL_cholesterol",
            #"feature_set": "LOFTEE_pLoF",
            #"feature_set": "AbExp_all_tissues",
            "feature_set": "LOFTEE_pLoF",
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

# %%
#Model with AGE+SEX+PC+PRS
restricted_variables = full_variables[1:34]

# %%
#Model with AGE+SEX+PC
basic_variables = restricted_variables[:-1]

# %%
model_variables = full_variables[35:]

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
#prs_features_df.write.parquet(snakemake.output["prs_features_pq"], mode="overwrite")

# %%
#prs_features_df = spark.read.parquet(snakemake.output["prs_features_pq"])

# %%
prs_features_pd_df = prs_features_df.toPandas()

# %%
spark.sparkContext._jvm.System.gc()

# %%
#Drop columns with no information
full_variables = prs_features_pd_df.columns[prs_features_pd_df.nunique() > 1]
prs_features_pd_df = prs_features_pd_df[full_variables]
prs_features_pd_df["sex_f31_0_0"] = prs_features_pd_df["sex_f31_0_0"].astype(int)

# %% [markdown]
# # read samples

# %%
train_test_split = pd.read_parquet(snakemake.input["train_test_split_pq"])
train_test_split

# %%
samples = train_test_split[["individual", "fold"]]
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
X_train = train_df[full_variables].drop(columns=[phenotype_col, "individual"])
y_train = train_df[phenotype_col]

# %%
X_test = test_df[full_variables].drop(columns=[phenotype_col,  "individual"])
y_test = test_df[phenotype_col]

# %% [markdown]
# ## train models

# %%
# Continuous traits:

# %%
### LGBMR

# %%
gbm_full = lgb.LGBMRegressor()
gbm_full.fit(X_train, y_train)
full_pred_test = gbm_full.predict(X_test)
full_model_r2 = sklearn.metrics.r2_score(y_test, full_pred_test)
full_model_r2

# %%
gbm_restricted = lgb.LGBMRegressor()
gbm_restricted.fit(X_train[restricted_variables], y_train)
restricted_pred_test = gbm_restricted.predict(X_test[restricted_variables])
restricted_model_r2 = sklearn.metrics.r2_score(y_test, restricted_pred_test)
restricted_model_r2

# %%
gbm_basic = lgb.LGBMRegressor()
gbm_basic.fit(X_train[basic_variables], y_train)
basic_pred_test = gbm_basic.predict(X_test[basic_variables])
basic_model_r2 = sklearn.metrics.r2_score(y_test, basic_pred_test)
basic_model_r2

# %%
# nr_of_percentiles = 10

pred_df = (
    test_df[["individual", "age_when_attended_assessment_centre_f21003_0_0", "sex_f31_0_0", phenotype_col]]
    .rename(columns={phenotype_col: "measurement", "age_when_attended_assessment_centre_f21003_0_0" : "age", "sex_f31_0_0": "sex"})
    .assign(**{
        "phenotype_col": phenotype_col,
        "method": snakemake.wildcards["feature_set"],
        # "phenotype_percentile": np.array(pd.qcut(test_df[phenotype_col], nr_of_quantiles, labels=False)),
        # "at_risk_high": np.array(pd.qcut(test_df[phenotype_col], 10, labels=False)) == nr_of_percentiles - 1,
        # "at_risk_low": np.array(pd.qcut(test_df[phenotype_col], 10, labels=False)) == 0,
        "full_model_pred": full_pred_test,
        # "full_model_pred_percentile": np.array(pd.qcut(full_pred_test, nr_of_quantiles, labels=False)),
        "restricted_model_pred": restricted_pred_test,
        # "restricted_model_pred_percentile": np.array(pd.qcut(restricted_pred_test, nr_of_quantiles, labels=False)),
        "basic_model_pred": basic_pred_test,
        #"basic_model_pred_quantile": np.array(pd.qcut(basic_pred_test, nr_of_quantiles, labels=False)),
        **snakemake.wildcards,
    })
)
pred_df

# %% [raw]
# pred_df = (
#     pred_df
#     .assign(**{
#         "residuals_full": pred_df["measurement"] - pred_df["full_model_pred"],
#         "residuals_restricted": pred_df["measurement"] - pred_df["restricted_model_pred"],
#         "residuals_basic": pred_df["measurement"] - pred_df["basic_model_pred"],
#         "abs_diff_full_restricted_pred": (pred_df["full_model_pred"] - pred_df["basic_model_pred"]).abs(),
#         "full_model_pred_rank": pred_df["full_model_pred"].rank(ascending = False),
#         "restricted_model_pred_rank": pred_df["restricted_model_pred"].rank(ascending = False),
#         "basic_model_pred_rank": pred_df["basic_model_pred"].rank(ascending = False),
#         #"pred_full_at_risk_low": pred_df["full_model_pred_quantile"] < 2,
#         #"pred_restricted_at_risk_low": pred_df["restricted_model_pred_quantile"] < 2
#     })
# )
# pred_df

# %%
# Save predictions
pred_df.to_parquet(snakemake.output["predictions_pq"], index=False)

# %% [raw]
# plot_df = pred_df[["phenotype_col", "measurement", "basic_model_pred", "restricted_model_pred", "full_model_pred"]].rename(columns={
#     "restricted_model_pred": f"Age+Sex+PC+PRS \n r²={restricted_model_r2:.3f}" ,"full_model_pred": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']} \n r²={full_model_r2:.3f}", "basic_model_pred": f"Age+Sex+PC \n r²={basic_model_r2:.3f}"
# })
# plot_df = plot_df.melt(id_vars=["phenotype_col", "measurement"])
#
# plot = (
#     pn.ggplot(plot_df, pn.aes(y="measurement", x="value"))
#     + pn.ggtitle(f"Predictions of LGBM models for {phenotype_col} ({phenocode})")
#     + pn.geom_bin_2d(bins=100)
#     + pn.geom_smooth(method="lm", color="red")
#     + pn.facet_grid("phenotype_col ~ variable")
#     + pn.scale_fill_continuous(trans = "log10")
#     + pn.theme(figure_size=(12, 4))
#     + pn.theme(title = pn.element_text(va = "top", linespacing = 4))
#     + pn.coord_equal()
#     + pn.labs(
#         x="prediction",
#     )
# )
#
# pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
# display(plot)

# %%
# Calc PRC
prc_df = []
for model in ["basic", "restricted", "full"]:
    for extreme in ["bottom", "top"]:
        for percentile in [0.01 , 0.05, 0.1, 0.2]:
            if extreme == "top":
                measurement_quantile = pred_df["measurement"].quantile(1-percentile)
                prc = sklearn.metrics.precision_recall_curve(pred_df["measurement"].ge(measurement_quantile), pred_df[f"{model}_model_pred"])
            else:
                measurement_quantile = pred_df["measurement"].quantile(percentile)
                prc = sklearn.metrics.precision_recall_curve(~pred_df["measurement"].ge(measurement_quantile), (-1) * pred_df[f"{model}_model_pred"])
            prc_df.append(pd.DataFrame({"extreme": extreme, "percentile": percentile, "method": model, "precision" : prc[0], "recall" : prc[1], "auPRC": sklearn.metrics.auc(prc[1], prc[0])}))
prc_df = pd.concat(prc_df)
prc_df["method"] = prc_df["method"].replace({"basic": "Age+Sex+PC", "restricted": "Age+Sex+PC+PRS", "full": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}"})

# %%
# Save prc 
prc_df.query("method=='Age+Sex+PC' or method=='Age+Sex+PC+PRS'").to_parquet(snakemake.output["precision_recall_baseline_pq"], index=False)
prc_df.query(f"method=='Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}'").to_parquet(snakemake.output["precision_recall_full_pq"], index=False)

# %%
prc_df

# %% [raw]
# plot = (
#     pn.ggplot(prc_df, pn.aes(x="recall", y="precision", color="method"))
#     + pn.geom_step()
#     + pn.facet_grid("extreme ~ percentile", labeller = 'label_both')
#     + pn.theme(figure_size=(16, 8))
#     + pn.ggtitle(f"Precision-Recall curves of LGBMRegressor models predicting extreme {phenotype_col} ({phenocode})")
# )
# pn.ggsave(plot = plot, filename = snakemake.output["prc_plot_png"], dpi=DPI)
# display(plot)

# %% [raw]
# risk_type = "at_risk_low"
# age_data_df = pd.DataFrame({
#     "baseline": pred_df.groupby("age")[risk_type].mean(), 
#     "top_decile_AbExp" : pred_df.query("full_model_pred_quantile==9").groupby("age")[risk_type].mean(),
#     "bottom_decile_AbExp" : pred_df.query("full_model_pred_quantile==0").groupby("age")[risk_type].mean(),
#     "top_decile_PRS" : pred_df.query("restricted_model_pred_quantile==9").groupby("age")[risk_type].mean(),
#     "bottom_decile_PRS" : pred_df.query("restricted_model_pred_quantile==0").groupby("age")[risk_type].mean()
# }).dropna()

# %% [raw]
# plot = (pn.ggplot(age_data_df.reset_index(), pn.aes(x='age')) +\
#     pn.geom_line(pn.aes(y='baseline'), color='blue') +\
#     pn.geom_line(pn.aes(y='bottom_decile_AbExp'), color='green') +\
#     pn.geom_line(pn.aes(y='top_decile_AbExp'), color='green') +\
#     pn.geom_line(pn.aes(y='bottom_decile_PRS'), color='red') +\
#     pn.geom_line(pn.aes(y='top_decile_PRS'), color='red')
# )
# display(plot)

# %% [raw] {"tags": [], "jp-MarkdownHeadingCollapsed": true}
# ## Testing

# %% [raw]
# ax = lgb.plot_importance(gbm_full, height=.5, max_num_features=40, title=f"Feature importance in LGBM model predicting {phenotype_col} (top 40 features)")
# ax.figure.savefig(f"../../plots/{phenotype_col}/LGBM_feature_importance.png")

# %% [raw]
# ax = pd.DataFrame({"sex+age+PC": pred_df.groupby("basic_model_pred_quantile")["phenotype_quantile"].mean(),"+PRS": pred_df.groupby("restricted_model_pred_quantile")["phenotype_quantile"].mean(), "+AbExpAllTissuesForSignificantGenes": pred_df.groupby("full_model_pred_quantile")["phenotype_quantile"].mean()}).plot()
# ax.set_ylabel("Measurement mean decile")
# ax.set_xlabel("Prediction decile")
# ax.figure.savefig(f"../../plots/{phenotype_col}/LGBM_mean_measurement_decile_by_prediction_decile.png")

# %% [raw]
# ax = pd.DataFrame({"sex+age+PC": pred_df.groupby("phenotype_quantile").apply(lambda decile: calculate_r2(decile["measurement"],decile["basic_model_pred"])),"+PRS": pred_df.groupby("phenotype_quantile").apply(lambda decile: calculate_r2(decile["measurement"],decile["restricted_model_pred"])),"+PRS+AbExpAllTissuesForSignificantGenes": pred_df.groupby("phenotype_quantile").apply(lambda decile: calculate_r2(decile["measurement"],decile["full_model_pred"])), "+PRS+Loftee": pred_df.groupby("phenotype_quantile").apply(lambda decile: calculate_r2(decile["measurement"],decile["loftee_model_pred"]))}).plot()
# ax.set_xlabel(f"{phenotype_col} measurement decile")
# ax.set_ylabel("r^2")
# #ax.figure.savefig(f"../../plots/{phenotype_col}/LGBM_r2_by_measurement_decile.png")

# %% [raw]
# ax = pd.DataFrame({"sex+age+PC": pred_df.groupby("phenotype_quantile").apply(lambda decile: sklearn.metrics.mean_squared_error(decile["measurement"],decile["basic_model_pred"])),"+PRS": pred_df.groupby("phenotype_quantile").apply(lambda decile: sklearn.metrics.mean_squared_error(decile["measurement"],decile["restricted_model_pred"])),"+AbExpAllTissuesForSignificantGenes": pred_df.groupby("phenotype_quantile").apply(lambda decile: sklearn.metrics.mean_squared_error(decile["measurement"],decile["full_model_pred"]))}).plot()
# ax.set_xlabel(f"{phenotype_col} measurement decile")
# ax.set_ylabel("MSE")
# #ax.figure.savefig(f"../../plots/{phenotype_col}/LGBM_MSE_by_measurement_decile.png")

# %% [raw]
# ax = pd.DataFrame({"sex+age+PC": pred_df.groupby("phenotype_quantile").apply(lambda decile: (decile["phenotype_quantile"]==decile["basic_model_pred_quantile"]).mean()),"+PRS": pred_df.groupby("phenotype_quantile").apply(lambda decile: (decile["phenotype_quantile"]==decile["restricted_model_pred_quantile"]).mean()),"+AbExpAllTissuesForSignificantGenes": pred_df.groupby("phenotype_quantile").apply(lambda decile: (decile["phenotype_quantile"]==decile["full_model_pred_quantile"]).mean()), "+Loftee": pred_df.groupby("phenotype_quantile").apply(lambda decile: (decile["phenotype_quantile"]==decile["loftee_model_pred_quantile"]).mean())}).plot()
# ax.set_xlabel("measurement decile")
# ax.set_ylabel("prediction decile match")
# ax.figure.savefig(f"../../plots/{phenotype_col}/LGBM_prediction_decile_accuracy_by_measurement_decile.png")

# %% [raw]
# models = ["loftee", "full", "restricted", "basic"]
# odd_ratios = {model: {} for model in models}
# for model in models:
#     for q in np.logspace(np.log10(50), np.log10(10000), num=50):
#         top_samples = pred_df.query(f"{model}_model_pred_rank<={q}")
#         bottom_samples = pred_df.query(f"{model}_model_pred_rank>{q}")
#         
#         #Version1:
#         #odd_ratios[model][q] = (top_samples["at_risk"].mean()) / ((bottom_samples["at_risk"]).mean())
#         
#         #Version2:
#         TP = top_samples["at_risk"].sum()
#         TN = (~bottom_samples["at_risk"]).sum()
#         FP = (~top_samples["at_risk"]).sum()
#         FN = bottom_samples["at_risk"].sum()
#         odd_ratios[model][q] = (TP * TN) / (FP * FN)
#         
#         #Both versions are equivalent
# odd_ratios = pd.DataFrame.from_dict(odd_ratios, orient="index").transpose()

# %% [raw]
# fig1, ax1 = plt.subplots()
# ax1.plot(odd_ratios["basic"], label="Age+Sex+PC")
# ax1.plot(odd_ratios["restricted"], label="+PRS")
# ax1.plot(odd_ratios["full"], label="+AbExpSignificantGenesAllTissues")
# ax1.plot(odd_ratios["loftee"], label="+Loftee")
# plt.legend()
# ax1.set_xscale('log')
# ax1.invert_xaxis()
# plt.xlabel("Top n Samples by prediction")
# plt.ylabel("Odds ratio")
# plt.title(f"Models predicting {phenotype_col} - Test for measurement in top decile")
# ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# ax1.figure.savefig(f"../../plots/{phenotype_col}/LGBM_OddsRatio_Top_Decile_Test.png")

# %% [raw]
# model_discrepency_df = pred_df.query("abs_diff_full_restricted_pred>1")
# sklearn.metrics.mean_squared_error(model_discrepency_df["measurement"], model_discrepency_df["full_model_pred"])

# %% [raw]
# sklearn.metrics.mean_squared_error(model_discrepency_df["measurement"], model_discrepency_df["restricted_model_pred"])

# %% [raw]
# # Binary traits:

# %% [raw]
# y_train

# %% [raw]
# model_pipeline = sklearn.pipeline.Pipeline([
#     ("normalize", sklearn.preprocessing.StandardScaler()),
#     ("logistic_regression_model", sklearn.linear_model.LogisticRegressionCV(n_jobs=16, cv=cv_split, max_iter=250, verbose=False, penalty = 'elasticnet', solver = 'saga', l1_ratios = [0.5])),
# ])

# %% [raw]
# full_model = model_pipeline.fit(X_train, y_train)
# full_model

# %% [raw]
# full_pred_test = full_model.predict(X_test)
# full_pred_test

# %% [raw]
# full_pred_prob_test = full_model.predict_proba(X_test)[::,1]
# full_pred_prob_test

# %% [raw]
# full_model.score(X_test, y_test)

# %% [raw]
# cm = sklearn.metrics.confusion_matrix(y_test, full_pred_test)

# %% [raw]
# cm

# %% [raw]
# model_pipeline_res = sklearn.pipeline.Pipeline([
#     ("normalize", sklearn.preprocessing.StandardScaler()),
#     ("logistic_regression_model", sklearn.linear_model.LogisticRegressionCV(n_jobs=16, cv=cv_split, max_iter=250, verbose=False, penalty = 'elasticnet', solver = 'saga', l1_ratios = [0.5])),
# ])

# %% [raw]
# restricted_model = model_pipeline_res.fit(X_train[restricted_variables], y_train)
# restricted_model

# %% [raw]
# restricted_pred_test = restricted_model.predict(X_test[restricted_variables])
# restricted_pred_test

# %% [raw]
# restricted_pred_prob_test = restricted_model.predict_proba(X_test[restricted_variables])[::,1]
# restricted_pred_prob_test

# %% [raw]
# restricted_model.score(X_test[restricted_variables], y_test)

# %% [raw]
# cm = sklearn.metrics.confusion_matrix(y_test, restricted_pred_test)

# %% [raw]
# cm

# %% [raw]
# model_pipeline_basic = sklearn.pipeline.Pipeline([
#     ("normalize", sklearn.preprocessing.StandardScaler()),
#     ("logistic_regression_model", sklearn.linear_model.LogisticRegressionCV(n_jobs=8, cv=cv_split, class_weight="balanced", max_iter=250, verbose=False, penalty = 'elasticnet', solver = 'saga', l1_ratios = [0.5])),
# ])

# %% [raw]
# basic_model = model_pipeline_basic.fit(X_train[basic_variables], y_train)
# basic_model

# %% [raw]
# basic_pred_test = basic_model.predict(X_test[basic_variables])
# basic_pred_test

# %% [raw]
# basic_pred_prob_test = basic_model.predict_proba(X_test[basic_variables])[::,1]
# basic_pred_prob_test

# %% [raw]
# basic_model.score(X_test[basic_variables], y_test)

# %% [raw]
# sklearn.metrics.confusion_matrix(y_test, basic_pred_test)

# %% [raw]
# y_train_binary = pd.qcut(y_train, 100, labels=False)>=5
# y_test_binary = pd.qcut(y_test, 100, labels=False)>=5

# %% [raw]
# gbm_classifier_full = lgb.LGBMClassifier()
# gbm_classifier_full.fit(X_train[[c for c in X_train.columns if "@pLof" not in c]], y_train_binary)
# full_pred_prob_test = gbm_classifier_full.predict_proba(X_test[[c for c in X_train.columns if "@pLof" not in c]])[::,1]
# gbm_classifier_full.score(X_test[[c for c in X_train.columns if "@pLof" not in c]], y_test_binary)

# %% [raw]
# gbm_classifier_loftee = lgb.LGBMClassifier()
# gbm_classifier_loftee.fit(X_train[[c for c in X_train.columns if "@AbExp" not in c]], y_train_binary)
# loftee_pred_prob_test = gbm_classifier_loftee.predict_proba(X_test[[c for c in X_train.columns if "@AbExp" not in c]])[::,1]
# gbm_classifier_loftee.score(X_test[[c for c in X_train.columns if "@AbExp" not in c]], y_test_binary)

# %% [raw]
# gbm_classifier_restricted = lgb.LGBMClassifier()
# gbm_classifier_restricted.fit(X_train[restricted_variables], y_train_binary)
# restricted_pred_prob_test = gbm_classifier_restricted.predict_proba(X_test[restricted_variables])[::,1]
# gbm_classifier_restricted.score(X_test[restricted_variables], y_test_binary)

# %% [raw]
# gbm_classifier_basic = lgb.LGBMClassifier()
# gbm_classifier_basic.fit(X_train[basic_variables], y_train_binary)
# basic_pred_prob_test = gbm_classifier_basic.predict_proba(X_test[basic_variables])[::,1]
# gbm_classifier_basic.score(X_test[basic_variables], y_test_binary)

# %% [raw]
# fpr_full, tpr_full, _ = sklearn.metrics.roc_curve(y_test_binary,  full_pred_prob_test)
# auc_full = sklearn.metrics.roc_auc_score(y_test_binary, full_pred_prob_test)
#
# fpr_loftee, tpr_loftee, _ = sklearn.metrics.roc_curve(y_test_binary,  loftee_pred_prob_test)
# auc_loftee = sklearn.metrics.roc_auc_score(y_test_binary, loftee_pred_prob_test)
#
# fpr_res, tpr_res, _ = sklearn.metrics.roc_curve(y_test_binary,  restricted_pred_prob_test)
# auc_res = sklearn.metrics.roc_auc_score(y_test_binary, restricted_pred_prob_test)
#
# fpr_basic, tpr_basic, _ = sklearn.metrics.roc_curve(y_test_binary,  basic_pred_prob_test)
# auc_basic = sklearn.metrics.roc_auc_score(y_test_binary, basic_pred_prob_test)
#
# #create ROC curve
# plt.plot(fpr_basic,tpr_basic,label="Age+Sex+PC: AUC="+str(auc_basic))
# plt.plot(fpr_res,tpr_res,label="+PRS: AUC="+str(auc_res))
# plt.plot(fpr_loftee,tpr_loftee,label="+PRS+Loftee: AUC="+str(auc_loftee))
# plt.plot(fpr_full,tpr_full,label="+PRS+AbExp: AUC="+str(auc_full))
# plt.ylabel('True Positive Rate')
# plt.xlabel('False Positive Rate')
# plt.legend(loc=4)
# plt.show()

# %% [raw]
# fig, ax = plt.subplots()
# sklearn.metrics.PrecisionRecallDisplay.from_predictions(y_test_binary, pred_df["basic_model_pred"], name="Age+Sex+PC:" ,ax=ax)
# sklearn.metrics.PrecisionRecallDisplay.from_predictions(y_test_binary, pred_df["restricted_model_pred"], name="+PRS:", ax=ax)
# sklearn.metrics.PrecisionRecallDisplay.from_predictions(y_test_binary, pred_df["loftee_model_pred"], name="+PRS+Loftee:",ax=ax)
# sklearn.metrics.PrecisionRecallDisplay.from_predictions(y_test_binary, pred_df["full_model_pred"], name="+PRS+AbExpAllTissuesForSignificantGenes:",ax=ax)
#
# plt.title(f"Precision-Recall curves of LGBMRegressor models predicting extreme {phenotype_col} ({phenocode})")
# plt.xlabel("Recall")
# plt.ylabel("Precision")
# plt.legend(loc=1)
# plt.show()

# %% [raw]
# # save model + predictions

# %% [raw]
# snakemake.output

# %% [raw]
# ## predictions

# %% [raw]
# pred_df.to_parquet(snakemake.output["predictions_pq"], index=False)

# %% [raw]
# ## summary

# %% [raw]
# restricted_model_summary = restricted_model.summary2()
# full_model_summary = full_model.summary2()

# %% [raw]
# with open(snakemake.output["restricted_summary_txt"], "w") as fd:
#     fd.write(restricted_model_summary.as_text())
#     fd.write("\n")

# %% [raw]
# with open(snakemake.output["full_summary_txt"], "w") as fd:
#     fd.write(full_model_summary.as_text())
#     fd.write("\n")

# %% [raw]
# ## params

# %% [raw]
# restricted_params_df = (
#     restricted_model_summary.tables[1]
#     .rename_axis(index="term")
#     .reset_index()
# )
# restricted_params_df

# %% [raw]
# restricted_params_df.to_parquet(snakemake.output["restricted_params_pq"], index=False)

# %% [raw]
# full_params_df = (
#     full_model_summary.tables[1]
#     .rename_axis(index="term")
#     .reset_index()
# )
# full_params_df

# %% [raw]
# full_params_df.to_parquet(snakemake.output["full_params_pq"], index=False)

# %% [raw]
# # plotting

# %% [raw]
# # TODO for Jonas

# %% [raw]
# ## phenotype correlation

# %% [raw]
# plot_df = pred_df[["phenotype_col", "measurement", "restricted_model_pred", "full_model_pred", "loftee_model_pred"]].rename(columns={
#     "restricted_model_pred": f"PRS r²={restricted_model_r2:.3f}" ,"full_model_pred": f"PRS+AbExp r²={full_model_r2:.3f}", "loftee_model_pred": f"PRS+Loftee r²={loftee_model_r2:.3f}"
# })
# plot_df = plot_df.melt(id_vars=["phenotype_col", "measurement"])
#
# plot = (
#     pn.ggplot(plot_df, pn.aes(y="measurement", x="value"))
#     + pn.ggtitle(f"Predictions for {phenotype_col} of models on PRS, PRS + AbExp and PRS + Loftee")
#     + pn.geom_bin_2d(bins=100)
#     + pn.geom_smooth(method="lm", color="red")
#     + pn.facet_grid("phenotype_col ~ variable")
#     + pn.scale_fill_continuous(trans = "log10")
#     + pn.coord_equal()
#     + pn.labs(
#         x="prediction",
#     )
# )
#
# pn.ggsave(plot = plot, filename = f"../../plots/{phenotype_col}/LGBM_Scatter_PRS_and_PRS+AbExp.png", dpi=DPI)
# display(plot)

# %% [raw]
# plot_df = pred_df[["full_model_pred", "restricted_model_pred"]].rename(columns = {"full_model_pred": "+AbExp (prediction)", "restricted_model_pred": "PRS (prediction)"})
#
# plot = (
#     pn.ggplot(plot_df, pn.aes(y="+AbExp (prediction)", x="PRS (prediction)"))
#     + pn.ggtitle(f"Predictions for {phenotype_col} of models on PRS vs. PRS + AbExp")
#     + pn.geom_bin_2d(bins=100)
#     + pn.geom_smooth(method="lm", color="red")
#     #+ pn.facet_grid("phenotype_col ~ variable")
#     + pn.scale_fill_continuous(trans = "log10")
#     #+ pn.coord_equal()
#     #+ pn.labs(
#     #    x="prediction",
#     #)
# )
# pn.ggsave(plot = plot, filename = f"../../plots/{phenotype_col}/LGBM_Scatter_PRS_vs_PRS+AbExp.png", dpi=DPI)
# display(plot)

# %% [raw]
# plot_df = pred_df[["loftee_model_pred", "restricted_model_pred"]].rename(columns = {"loftee_model_pred": "+loftee (prediction)", "restricted_model_pred": "PRS (prediction)"})
#
# plot = (
#     pn.ggplot(plot_df, pn.aes(y="+loftee (prediction)", x="PRS (prediction)"))
#     + pn.ggtitle(f"Predictions for {phenotype_col} of models on PRS vs. PRS + Loftee")
#     + pn.geom_bin_2d(bins=100)
#     + pn.geom_smooth(method="lm", color="red")
#     #+ pn.facet_grid("phenotype_col ~ variable")
#     + pn.scale_fill_continuous(trans = "log10")
#     #+ pn.coord_equal()
#     #+ pn.labs(
#     #    x="prediction",
#     #)
# )
# pn.ggsave(plot = plot, filename = f"../../plots/{phenotype_col}/LGBM_Scatter_PRS_vs_PRS+Loftee.png", dpi=DPI)
# display(plot)

# %% [raw]
# plot_df = pred_df[["loftee_model_pred", "full_model_pred"]].rename(columns = {"loftee_model_pred": "PRS+loftee (prediction)", "full_model_pred": "PRS+AbExp (prediction)"})
#
# plot = (
#     pn.ggplot(plot_df, pn.aes(y="PRS+loftee (prediction)", x="PRS+AbExp (prediction)"))
#     + pn.ggtitle(f"Predictions for {phenotype_col} of models on PRS + AbExp vs. PRS + Loftee")
#     + pn.geom_bin_2d(bins=100)
#     + pn.geom_smooth(method="lm", color="red")
#     #+ pn.facet_grid("phenotype_col ~ variable")
#     + pn.scale_fill_continuous(trans = "log10")
#     #+ pn.coord_equal()
#     #+ pn.labs(
#     #    x="prediction",
#     #)
# )
# pn.ggsave(plot = plot, filename = f"../../plots/{phenotype_col}/LGBM_Scatter_PRS+AbExp_vs_PRS+Loftee.png", dpi=DPI)
# display(plot)

# %% [raw]
# test_df

# %% [raw]
# abexp_variables = [c for c in test_df.columns if "@AbExp" in c]
# plot_df = test_df[[c for c in test_df.columns if c in ["individual"]+abexp_variables]]
# plot_df = pd.melt(plot_df, id_vars=["individual"], value_vars=list(set(abexp_variables).intersection(set(plot_df.columns)))).rename(columns={"value": "score"})
# plot_df[["gene_id", "score_type"]] = plot_df["variable"].str.split("@",expand=True)
# plot_df["gene_symbol"] = plot_df["gene_id"]
# plot_df = pd.merge(plot_df, pred_df[["individual", "residuals_restricted"]])

# %% [raw]
# plot = (
#     pn.ggplot(plot_df, pn.aes(x="score", y="residuals_restricted"))
#     + pn.stat_bin2d(bins=100)
#     + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
#     + pn.facet_grid("gene_symbol ~ score_type", scales="free_x")
# #     + pn.labs(
# #         title=f"{phenotype_col}"
# #     )
#     + pn.theme(aspect_ratio=1/2, figure_size = (80, 160))
# )
# plot

# %% [raw]
# plot_max_df = plot_df.loc[plot_df.groupby(["individual", "gene_id"])['score'].idxmax()].rename(columns={"residuals_restricted": "residual"})
# plot_max_df["score_type"] = "AbExp_max"
# plot = (
#     pn.ggplot(plot_max_df, pn.aes(x="score", y="residual"))
#     + pn.stat_bin2d(bins=100)
#     + pn.geom_smooth(method = "lm", color="red")#, se = FALSE)
#     + pn.facet_grid("gene_symbol ~ score_type", scales="free_x")
#     + pn.theme(aspect_ratio=1/2)
#     + pn.ggtitle(f"Residuals of common PRS and {phenotype_col} vs. AbExpMax")
#     + pn.theme(title = pn.element_text(va = "top", linespacing = 8))
# )
# pn.ggsave(plot = plot, filename = f"../../plots/{phenotype_col}/LGBM_PRS_residuals_vs_AbExpMax.png", dpi=DPI)
# plot

# %% [raw]
# snakemake.params["output_basedir"]

# %% [raw]
# path = snakemake.params["output_basedir"] + "/phenotype_correlation"
# pn.ggsave(plot, path + ".png", dpi=DPI)
# pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [raw]
# from IPython.display import Image
# display(Image(snakemake.params["output_basedir"] + "/phenotype_correlation.png"))

# %% [raw]
# ## residual correlation

# %% [raw]
# snakemake.wildcards

# %% [raw]
# import sklearn.metrics as skmetrics
# import textwrap
#
# plot_r2 = skmetrics.r2_score(pred_df["residuals"], pred_df["full - restricted"])
#
# plot = (
#     pn.ggplot(pred_df, pn.aes(y="residuals", x="full - restricted"))
#     + pn.geom_bin_2d(bins=100)
#     + pn.geom_smooth(method="lm", color="red")
#     + pn.scale_fill_continuous(trans = "log10")
#     + pn.labs(
#         x="prediction\n(full model - restricted model)",
#         y=f"residuals\n(phenotype - restricted model)",
#         title=textwrap.dedent(
#             f"""
#             Polygenic score predicting '{phenotype_col}'
#             Covariates: '{snakemake.wildcards['covariates'].replace('_', '+')}'
#             r2: '{"%.4f" % plot_r2}'
#             (Every point represents one individual)\
#             """
#         ).strip(),
#         fill='Nr. of individuals'
#     )
#     + pn.theme(title=pn.element_text(linespacing=1.4))
#     + pn.coord_equal()
# )

# %% [raw]
# if len(restricted_variables) != len(full_variables):
#     display(plot)

# %% [raw]
# snakemake.params["output_basedir"]

# %% [raw]
# path = snakemake.params["output_basedir"] + "/residual_correlation"
# if len(restricted_variables) != len(full_variables):
#     pn.ggsave(plot, path + ".png", dpi=DPI)
#     pn.ggsave(plot, path + ".pdf", dpi=DPI)

# %% [raw]
#
