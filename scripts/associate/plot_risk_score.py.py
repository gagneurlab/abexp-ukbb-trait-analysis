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

import sklearn.metrics

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

from rep.notebook_init import setup_plot_style
setup_plot_style()

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %%
del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__polygenic_risk_score',
        default_wildcards={
            #"phenotype_col": "standing_height",
            #"phenotype_col": "glycated_haemoglobin_hba1c",
            #"phenotype_col": "Lipoprotein_A",
            #"phenotype_col": "BodyMassIndex",
            #"phenotype_col": "Triglycerides",
            "phenotype_col": "LDL_direct",
            #"phenotype_col": "systolic_blood_pressure",
            #"phenotype_col": "HDL_cholesterol",
            #"feature_set": "LOFTEE_pLoF",
            "feature_set": "AbExp_all_tissues",
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

# %% [markdown]
# ## Read Predictions

# %%
pred_df = pd.read_parquet(snakemake.output["predictions_pq"])
pred_df

# %%
full_model_r2 = sklearn.metrics.r2_score(pred_df["measurement"], pred_df["full_model_pred"])
restricted_model_r2 = sklearn.metrics.r2_score(pred_df["measurement"], pred_df["restricted_model_pred"])
basic_model_r2 = sklearn.metrics.r2_score(pred_df["measurement"], pred_df["basic_model_pred"])

# %%
nr_of_quantiles = 100

pred_df = (
    pred_df
    .assign(**{
        "phenotype_col": phenotype_col,
        "method": snakemake.wildcards["feature_set"],
        "phenotype_quantile": np.array(pd.qcut(pred_df["measurement"], nr_of_quantiles, labels=False)),
        "full_model_pred_quantile": np.array(pd.qcut(pred_df["full_model_pred"], nr_of_quantiles, labels=False)),
        "restricted_model_pred_quantile": np.array(pd.qcut(pred_df["restricted_model_pred"], nr_of_quantiles, labels=False)),
        "basic_model_pred_quantile": np.array(pd.qcut(pred_df["basic_model_pred"], nr_of_quantiles, labels=False)),
    })
)

# %%
pred_df = (
    pred_df
    .assign(**{
        "residuals_full": pred_df["measurement"] - pred_df["full_model_pred"],
        "residuals_restricted": pred_df["measurement"] - pred_df["restricted_model_pred"],
        "residuals_basic": pred_df["measurement"] - pred_df["basic_model_pred"],
        "abs_diff_full_restricted_pred": (pred_df["full_model_pred"] - pred_df["basic_model_pred"]).abs(),
        "full_model_pred_rank": pred_df["full_model_pred"].rank(ascending = False),
        "restricted_model_pred_rank": pred_df["restricted_model_pred"].rank(ascending = False),
        "basic_model_pred_rank": pred_df["basic_model_pred"].rank(ascending = False),
        "at_risk_low": pred_df["phenotype_quantile"] == 0,
        "at_risk_high": pred_df["phenotype_quantile"] == nr_of_quantiles-1
    })
)
pred_df

# %% [markdown]
# ## Read PRC

# %%
# Read prc 
prc_baseline_df = pd.read_parquet(snakemake.output["precision_recall_baseline_pq"])
prc_full_df = pd.read_parquet(snakemake.output["precision_recall_full_pq"])
prc_plof_df = pd.read_parquet("/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/LDL_direct/cov=sex_age_genPC_CLMP_PRS/fset=LOFTEE_pLoF/polygenic_risk_score/precision_recall.full.parquet")
prc_df = pd.concat([prc_baseline_df, prc_full_df, prc_plof_df])

# %% [markdown] {"tags": []}
# ## Plots

# %% [markdown] {"tags": []}
# ### Scatter Predictions

# %%
plot_df = pred_df[["phenotype_col", "measurement", "basic_model_pred", "restricted_model_pred", "full_model_pred"]].rename(columns={
    "restricted_model_pred": f"Age+Sex+PC+PRS \n r²={restricted_model_r2:.3f}" ,"full_model_pred": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']} \n r²={full_model_r2:.3f}", "basic_model_pred": f"Age+Sex+PC \n r²={basic_model_r2:.3f}"
})
plot_df = plot_df.melt(id_vars=["phenotype_col", "measurement"])

plot = (
    pn.ggplot(plot_df, pn.aes(y="measurement", x="value"))
    + pn.ggtitle(f"Predictions of LGBM models for {phenotype_col} ({phenocode})")
    + pn.geom_bin_2d(bins=100)
    + pn.geom_smooth(method="lm", color="red")
    + pn.facet_grid("phenotype_col ~ variable")
    + pn.scale_fill_continuous(trans = "log10")
    + pn.theme(figure_size=(12, 4))
    + pn.theme(title = pn.element_text(va = "top", linespacing = 4))
    + pn.coord_equal()
    + pn.labs(
        x="prediction",
    )
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
pred_df.query("full_model_pred_quantile==0")["at_risk_low"].agg(["size", "sum"])

# %%
pred_df.query("restricted_model_pred_quantile==0")["at_risk_low"].agg(["size", "sum"])

# %%
plot_df = pred_df[["full_model_pred", "restricted_model_pred"]].rename(columns = {"full_model_pred": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}", "restricted_model_pred": "Age+Sex+PC+PRS"})

plot = (
    pn.ggplot(plot_df, pn.aes(y=f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}", x="Age+Sex+PC+PRS"))
    + pn.ggtitle(f"Predictions for {phenotype_col} of models on PRS vs. PRS + {snakemake.wildcards['feature_set']}")
    + pn.geom_bin_2d(bins=100)
    + pn.geom_smooth(method="lm", color="red")
    + pn.scale_fill_continuous(trans = "log10")
    #+ pn.coord_equal()
)
display(plot)

# %% [markdown]
# ### Plot PRC for extremes

# %%
plot = (
    pn.ggplot(prc_df, pn.aes(x="recall", y="precision", color="method"))
    + pn.geom_step()
    + pn.facet_grid("percentile ~ extreme", labeller = 'label_both')
    + pn.theme(figure_size=(8, 8))
    + pn.ggtitle(f"Precision-Recall curves of LGBMRegressor models predicting extreme {phenotype_col} ({phenocode})")
    + pn.coord_equal()
)
#pn.ggsave(plot = plot, filename = snakemake.output["prc_plot_png"], dpi=DPI)
display(plot)

# %%
df = prc_df.groupby(["extreme", "percentile", "method"])["auPRC"].first()

# %%
df.reset_index()

# %%
plot = (
    pn.ggplot(df.reset_index(), pn.aes(y="auPRC", x="method"))
    + pn.ggtitle(f"Predictions of LGBM models for {phenotype_col} ({phenocode})")
    + pn.geom_bar(stat="identity")
    + pn.facet_grid("percentile ~ extreme")
    + pn.theme(figure_size=(6, 8))
    #+ pn.theme(axis_text_x = pn.element_text(rotation=90))
    + pn.coord_flip()
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
risk_type = "at_risk_low"
age_data_df = pd.DataFrame({
    "baseline": pred_df.groupby("age")[risk_type].mean(), 
    f"top_decile_PRS+{snakemake.wildcards['feature_set']}" : pred_df.query("full_model_pred_quantile==9").groupby("age")[risk_type].mean(),
    f"bottom_decile_PRS+{snakemake.wildcards['feature_set']}" : pred_df.query("full_model_pred_quantile==0").groupby("age")[risk_type].mean(),
    "top_decile_PRS" : pred_df.query("restricted_model_pred_quantile==9").groupby("age")[risk_type].mean(),
    "bottom_decile_PRS" : pred_df.query("restricted_model_pred_quantile==0").groupby("age")[risk_type].mean()
}).dropna().reset_index()
age_data_df = pd.melt(age_data_df, id_vars=['age'], value_vars=['baseline', f"top_decile_PRS+{snakemake.wildcards['feature_set']}", f"bottom_decile_PRS+{snakemake.wildcards['feature_set']}", 'top_decile_PRS', 'bottom_decile_PRS'], var_name='method', value_name='prevalence')

age_data_df["age"] = age_data_df["age"]

plot = (pn.ggplot(age_data_df, pn.aes(x='age', y='prevalence', color='method'))
    + pn.geom_line()
    + pn.ggtitle(f"Incidence of bottom decile {phenotype_col} ({phenocode}) for the top/bottom deciles of predictions")
)
display(plot)

# %%
risk_type = "at_risk_high"
age_data_df = pd.DataFrame({
    "baseline": pred_df.groupby("age")[risk_type].mean(), 
    f"top_decile_PRS+{snakemake.wildcards['feature_set']}" : pred_df.query("full_model_pred_quantile==9").groupby("age")[risk_type].mean(),
    f"bottom_decile_PRS+{snakemake.wildcards['feature_set']}" : pred_df.query("full_model_pred_quantile==0").groupby("age")[risk_type].mean(),
    "top_decile_PRS" : pred_df.query("restricted_model_pred_quantile==9").groupby("age")[risk_type].mean(),
    "bottom_decile_PRS" : pred_df.query("restricted_model_pred_quantile==0").groupby("age")[risk_type].mean()
}).dropna().reset_index()
age_data_df = pd.melt(age_data_df, id_vars=['age'], value_vars=['baseline', f"top_decile_PRS+{snakemake.wildcards['feature_set']}", f"bottom_decile_PRS+{snakemake.wildcards['feature_set']}", 'top_decile_PRS', 'bottom_decile_PRS'], var_name='method', value_name='prevalence')

plot = (pn.ggplot(age_data_df, pn.aes(x='age', y='prevalence', color='method'))
    + pn.geom_line()
    + pn.ggtitle(f"Incidence of top decile {phenotype_col} ({phenocode}) for the top/bottom deciles of predictions")
)
display(plot)

# %% [markdown]
# ### Plot odds-ratios (not very promising)

# %%
risk_type = "at_risk_high"
models = ["full", "restricted", "basic"]
odd_ratios = {model: {} for model in models}
for model in models:
    for q in np.logspace(np.log10(50), np.log10(10000), num=50):
        top_samples = pred_df.query(f"{model}_model_pred_rank<={q}")
        bottom_samples = pred_df.query(f"{model}_model_pred_rank>{q}")
        
        #Version1:
        #odd_ratios[model][q] = (top_samples["at_risk"].mean()) / ((bottom_samples["at_risk"]).mean())
        
        #Version2:
        TP = top_samples[risk_type].sum()
        TN = (~bottom_samples[risk_type]).sum()
        FP = (~top_samples[risk_type]).sum()
        FN = bottom_samples[risk_type].sum()
        odd_ratios[model][q] = (TP * TN) / (FP * FN)
        
        #Both versions are equivalent
odd_ratios = pd.DataFrame.from_dict(odd_ratios, orient="index").transpose()

# %%
fig1, ax1 = plt.subplots()
ax1.plot(odd_ratios["basic"], label="Age+Sex+PC")
ax1.plot(odd_ratios["restricted"], label="+PRS")
ax1.plot(odd_ratios["full"], label="+AbExpSignificantGenesAllTissues")
plt.legend()
ax1.set_xscale('log')
ax1.invert_xaxis()
plt.xlabel("Top n Samples by prediction")
plt.ylabel("Odds ratio")
plt.title(f"Models predicting {phenotype_col} - Test for measurement in top decile")
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
