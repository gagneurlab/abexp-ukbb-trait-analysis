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
            #"phenotype_col": "LDL_direct",
            #"phenotype_col": "systolic_blood_pressure",
            #"phenotype_col": "HDL_cholesterol",
            "phenotype_col": "Alanine_aminotransferase",
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

restricted_model_mean = pred_df["restricted_model_pred"].mean()
restricted_model_std = pred_df["restricted_model_pred"].std()

pehnotype_mean = pred_df["measurement"].mean()
phenotype_std = pred_df["measurement"].std()

# %%
nr_of_quantiles = 10

pred_df = (
    pred_df
    .assign(**{
        "phenotype_col": phenotype_col,
        "method": snakemake.wildcards["feature_set"],
        #"phenotype_quantile": np.array(pd.qcut(pred_df["measurement"], nr_of_quantiles, labels=False)),
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
        #"at_risk_low": pred_df["phenotype_quantile"] == 0,
        #"at_risk_high": pred_df["phenotype_quantile"] == nr_of_quantiles-1,
        #"at_risk": np.abs(pred_df["measurement"] - pehnotype_mean) > 1 * phenotype_std,
        "full_model_new_risk": np.abs(pred_df["full_model_pred"] - pred_df["restricted_model_pred"]) > (1 * restricted_model_std)
    })
)
pred_df

# %% [markdown]
# ## Read PRC

# %%
# Read prc 
prc_baseline_df = pd.read_parquet(snakemake.output["precision_recall_baseline_pq"])
prc_full_df = pd.read_parquet(snakemake.output["precision_recall_full_pq"])
prc_plof_df = pd.read_parquet("/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/Triglycerides/cov=sex_age_genPC_CLMP_PRS/fset=LOFTEE_pLoF/polygenic_risk_score/precision_recall.full.parquet")
prc_df = pd.concat([prc_baseline_df, prc_full_df])

# %% [markdown] {"tags": []}
# ## Plots

# %% [markdown] {"tags": []}
# ### Scatter Predictions

# %%
plot_df = pred_df[["phenotype_col", "measurement", "basic_model_pred", "restricted_model_pred", "full_model_pred", "full_model_new_risk"]].rename(columns={
    "restricted_model_pred": f"Age+Sex+PC+PRS \n r²={restricted_model_r2:.3f}" ,"full_model_pred": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']} \n r²={full_model_r2:.3f}", "basic_model_pred": f"Age+Sex+PC \n r²={basic_model_r2:.3f}"
})
plot_df = plot_df.melt(id_vars=["phenotype_col", "measurement", "full_model_new_risk"])

plot = (
    pn.ggplot(plot_df, pn.aes(y="measurement", x="value", color='full_model_new_risk'))
    + pn.ggtitle(f"Predictions of LGBM models for {phenotype_col} ({phenocode})")
    + pn.geom_bin_2d(bins=100)
    + pn.geom_smooth(method="lm", color="red")
    + pn.facet_grid("phenotype_col ~ variable")
    + pn.scale_fill_continuous(trans = "log10")
    + pn.theme(figure_size=(12, 4))
    + pn.theme(title = pn.element_text(va = "top", linespacing = 4))
    + pn.coord_equal()
    + pn.scale_color_manual(values = {True: "orange", False: "None"})
    + pn.labs(
        x="prediction",
    )
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
plot_df = pred_df[["full_model_pred", "restricted_model_pred", "full_model_new_risk"]].rename(columns = {"full_model_pred": f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}", "restricted_model_pred": "Age+Sex+PC+PRS"})

plot = (
    pn.ggplot(plot_df, pn.aes(y=f"Age+Sex+PC+PRS+{snakemake.wildcards['feature_set']}", x="Age+Sex+PC+PRS", color="full_model_new_risk"))
    + pn.ggtitle(f"Predictions for {phenotype_col} of models on PRS vs. PRS + {snakemake.wildcards['feature_set']}")
    + pn.geom_bin_2d(bins=100)
    + pn.geom_smooth(method="lm", color="red")
    + pn.scale_fill_continuous(trans = "log10")
    + pn.scale_color_manual(values = {True: "orange", False: "None"})
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
    + pn.theme(title = pn.element_text(va = "top", linespacing = 4))
    + pn.coord_equal()
)
#pn.ggsave(plot = plot, filename = snakemake.output["prc_plot_png"], dpi=DPI)
display(plot)

# %%
plot = (
    pn.ggplot(prc_df.groupby(["extreme", "percentile", "method"])["auPRC"].first().reset_index(), pn.aes(y="auPRC", x="method"))
    + pn.ggtitle(f"auPRC of LGBM models prediction extreme {phenotype_col} ({phenocode})")
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

# %%
res = []
for phenotype in ["LDL_direct", "standing_height", "glycated_haemoglobin_hba1c","Lipoprotein_A","BodyMassIndex","Triglycerides", "systolic_blood_pressure","HDL_cholesterol"]:
    pred_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=AbExp_all_tissues/polygenic_risk_score/predictions.parquet")
    pred_Loftee_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=LOFTEE_pLoF/polygenic_risk_score/predictions.parquet")
    pred_df["loftee_model_pred"] = pred_Loftee_df["full_model_pred"]
    
    restricted_model_mean = pred_df["restricted_model_pred"].mean()
    restricted_model_std = pred_df["restricted_model_pred"].std()

    pehnotype_mean = pred_df["measurement"].mean()
    phenotype_std = pred_df["measurement"].std()
    
    pred_df["at_risk"] = np.abs(pred_df["measurement"] - pehnotype_mean) > 2 * phenotype_std
    pred_df["full_model_new_risk"] = np.abs(pred_df["full_model_pred"] - restricted_model_mean) > 3 * restricted_model_std
    pred_df["loftee_model_new_risk"] = np.abs(pred_df["loftee_model_pred"] - restricted_model_mean) > 3 * restricted_model_std

    for method in ["full", "loftee"]:
        for metric in ["sum", "size", "mean"]: 
            res.append({
                "phenotype": phenotype,
                "method": method,
                "metric": metric,
                "individuals": pred_df.query(f"{method}_model_new_risk==True")["at_risk"].agg([metric])[metric]
            })
res_df = pd.DataFrame.from_dict(res)

# %%
res_df

# %%
plot = (
    pn.ggplot(res_df, pn.aes(y="individuals", x="method"))
    + pn.ggtitle(f"Additional number of at risk ({1/nr_of_quantiles} percentile) individuals identified compared to Age+Sex+PC+PRS")
    + pn.geom_bar(stat="identity")
    + pn.facet_grid("metric ~ phenotype", scales="free_y")
    #+ pn.facet_wrap(["phenotype", "extreme", "metric"])
    + pn.theme(figure_size=(16, 6))
    #+ pn.theme(axis_text_y = pn.element_text(rotation=90))
    #+ pn.coord_flip()
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
res = []
nr_of_quantiles = 20
for phenotype in ["LDL_direct", "standing_height", "glycated_haemoglobin_hba1c","Lipoprotein_A","BodyMassIndex","Triglycerides", "systolic_blood_pressure","HDL_cholesterol"]:
    pred_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=AbExp_all_tissues/polygenic_risk_score/predictions.parquet")
    pred_Loftee_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=LOFTEE_pLoF/polygenic_risk_score/predictions.parquet")
    pred_df["loftee_model_pred"] = pred_Loftee_df["full_model_pred"]
    pred_df["phenotype_quantile"] = pd.qcut(pred_df["measurement"], nr_of_quantiles, labels=False, duplicates="drop")
    pred_df["restricted_model_pred_quantile"] = pd.qcut(pred_df["restricted_model_pred"], nr_of_quantiles, labels=False, duplicates="drop")
    pred_df["full_model_pred_quantile"] = pd.qcut(pred_df["full_model_pred"], nr_of_quantiles, labels=False, duplicates="drop")
    pred_df["loftee_model_pred_quantile"] = pd.qcut(pred_df["loftee_model_pred"], nr_of_quantiles, labels=False, duplicates="drop")
    pred_df["at_risk_low"] = (pred_df["phenotype_quantile"] == 0)
    pred_df["at_risk_high"] = (pred_df["phenotype_quantile"] == nr_of_quantiles-1)
    for extreme in ["low", "high"]:
        for method in ["full", "loftee"]:
            for metric in ["sum", "size", "mean"]: 
                res.append({
                    "phenotype": phenotype,
                    "extreme": extreme,
                    "method": method,
                    "metric": metric,
                    "individuals": pred_df.query(f"{method}_model_pred_quantile=={0 if extreme == 'low' else nr_of_quantiles-1}")[f"at_risk_{extreme}"].agg([metric])[metric] - pred_df.query(f"restricted_model_pred_quantile=={0 if extreme == 'low' else nr_of_quantiles-1}")[f"at_risk_{extreme}"].agg([metric])[metric]
                })
res_df = pd.DataFrame.from_dict(res)

# %%
plot = (
    pn.ggplot(res_df.query("metric=='sum'"), pn.aes(y="individuals", x="method"))
    + pn.ggtitle(f"Additional number of at risk ({1/nr_of_quantiles} percentile) individuals identified compared to Age+Sex+PC+PRS")
    + pn.geom_bar(stat="identity")
    + pn.facet_grid("extreme ~ phenotype", scales="free_y")
    #+ pn.facet_wrap(["phenotype", "extreme", "metric"])
    + pn.theme(figure_size=(16, 6))
    #+ pn.theme(axis_text_y = pn.element_text(rotation=90))
    #+ pn.coord_flip()
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
res = []
nr_of_quantiles = 10
distance_std = 0.5
for phenotype in ["LDL_direct", "standing_height", "glycated_haemoglobin_hba1c","Lipoprotein_A","BodyMassIndex","Triglycerides","HDL_cholesterol"]:
    pred_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=AbExp_all_tissues/polygenic_risk_score/predictions.parquet")
    pred_Loftee_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=LOFTEE_pLoF/polygenic_risk_score/predictions.parquet")
    pred_df["loftee_model_pred"] = pred_Loftee_df["full_model_pred"]
    pred_df["abexp_model_pred"] = pred_df["full_model_pred"]
    pred_df["phenotype_quantile"] = pd.qcut(pred_df["measurement"], nr_of_quantiles, labels=False, duplicates="drop")
    pred_df["at_risk_low"] = (pred_df["phenotype_quantile"] == 0)
    pred_df["at_risk_high"] = (pred_df["phenotype_quantile"] == nr_of_quantiles-1)
    distance_to_common_prs = pred_df["restricted_model_pred"].std() * distance_std
    print(phenotype, distance_to_common_prs)
    
    for extreme in ["low", "high"]:
        for method in ["abexp", "loftee"]:
            pred_df[f"new_{method}_at_risk_{extreme}"] = ((pred_df["restricted_model_pred"] - pred_df[f"{method}_model_pred"]) > distance_to_common_prs) if extreme == 'low' else ((pred_df["restricted_model_pred"] - pred_df[f"{method}_model_pred"]) < (-distance_to_common_prs))
            for model_type in ["baseline", "full"]:
                res.append({
                    "phenotype": phenotype,
                    "extreme": extreme,
                    "method": method,
                    "model_type": model_type,
                    "individuals": pred_df[f"new_{method}_at_risk_{extreme}"].sum(),
                    "mse": sklearn.metrics.mean_squared_error(pred_df.query(f"new_{method}_at_risk_{extreme}")["measurement"], pred_df.query(f"new_{method}_at_risk_{extreme}")[f"{method}_model_pred" if model_type=="full" else "restricted_model_pred"]) if pred_df[f"new_{method}_at_risk_{extreme}"].sum()>0 else 0,
                    "rmse": sklearn.metrics.mean_squared_error(pred_df.query(f"new_{method}_at_risk_{extreme}")["measurement"], pred_df.query(f"new_{method}_at_risk_{extreme}")[f"{method}_model_pred" if model_type=="full" else "restricted_model_pred"], squared=False) if pred_df[f"new_{method}_at_risk_{extreme}"].sum()>0 else 0,
                    "r2": sklearn.metrics.r2_score(pred_df.query(f"new_{method}_at_risk_{extreme}")["measurement"], pred_df.query(f"new_{method}_at_risk_{extreme}")[f"{method}_model_pred" if model_type=="full" else "restricted_model_pred"]) if pred_df[f"new_{method}_at_risk_{extreme}"].sum()>0 else 0

                })
res_df = pd.DataFrame.from_dict(res)

# %%
plot = (
    pn.ggplot(res_df, pn.aes(y="individuals", x="method"))
    + pn.ggtitle(f"Number of individuals where prediction differs by more than {distance_std} standard deviation(s) from common PRS")
    + pn.geom_bar(stat="identity")
    + pn.facet_grid("extreme ~ phenotype")
    #+ pn.facet_wrap(["phenotype", "extreme", "metric"])
    + pn.theme(figure_size=(16, 6))
    #+ pn.theme(axis_text_y = pn.element_text(rotation=90))
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
plot = (
    pn.ggplot(res_df, pn.aes(y="rmse", x="method", fill="model_type"))
    + pn.ggtitle(f"RMSE of predictions for individuals where prediction differs by more than {distance_std} standard deviation(s) from common PRS")
    + pn.geom_bar(stat="identity", position="dodge")
    + pn.facet_grid("extreme ~ phenotype")
    #+ pn.facet_wrap(["phenotype", "extreme", "metric"])
    + pn.theme(figure_size=(16, 6))
    #+ pn.theme(axis_text_y = pn.element_text(rotation=90))
    #+ pn.coord_flip()
)

#pn.ggsave(plot = plot, filename = snakemake.output["predictions_plot_png"], dpi=DPI)
display(plot)

# %%
res = []
nr_of_quantiles = 100
distance_std = 1
for phenotype in ["LDL_direct", "standing_height", "glycated_haemoglobin_hba1c","Lipoprotein_A","BodyMassIndex","Triglycerides","HDL_cholesterol"]:
    pred_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=AbExp_all_tissues/polygenic_risk_score/predictions.parquet")
    pred_Loftee_df = pd.read_parquet(f"/s/project/rep/processed/trait_associations_v3/ukbb_wes_200k/associate/{phenotype}/cov=sex_age_genPC_CLMP_PRS/fset=LOFTEE_pLoF/polygenic_risk_score/predictions.parquet")
    pred_df["LOFTEE_model_pred"] = pred_Loftee_df["full_model_pred"]
    pred_df["AbExp_model_pred"] = pred_df["full_model_pred"]
    pred_df["phenotype_quantile"] = pd.qcut(pred_df["measurement"], nr_of_quantiles, labels=False, duplicates="drop")
    pred_df["at_risk_low"] = (pred_df["phenotype_quantile"] == 0)
    pred_df["at_risk_high"] = (pred_df["phenotype_quantile"] == nr_of_quantiles-1)
    pred_df["at_risk"] = pred_df["at_risk_high"] | pred_df["at_risk_low"]
    distance_to_common_prs = pred_df["restricted_model_pred"].std() * distance_std
    print(phenotype, distance_to_common_prs)
    
    for method in ["AbExp", "LOFTEE"]:
        pred_df[f"{method}_deviation_from_PRS"] = np.abs((pred_df["restricted_model_pred"] - pred_df[f"{method}_model_pred"]))
        pred_df[f"new_{method}_at_risk"] = (pred_df[f"{method}_deviation_from_PRS"]) > distance_to_common_prs
        #pred_df[f"new_{method}_at_risk"] = pd.qcut(pred_df[f"{method}_deviation_from_PRS"], nr_of_quantiles, labels=False, duplicates="drop").isin([0,nr_of_quantiles-1])
        for model_type in ["common_PRS", "full"]:
            for at_risk in [True, False]:
                res.append({
                    "phenotype": phenotype,
                    "method": method,
                    #f"extreme_phenotype_{1/nr_of_quantiles*100}%": at_risk,
                    "extreme_phenotype": at_risk,
                    "model_type": model_type,
                    "individuals": pred_df[f"new_{method}_at_risk"].sum(),
                    "mse": sklearn.metrics.mean_squared_error(pred_df.query(f"new_{method}_at_risk")["measurement"], pred_df.query(f"new_{method}_at_risk")[f"{method}_model_pred" if model_type=="full" else "restricted_model_pred"]) if pred_df[f"new_{method}_at_risk"].sum()>0 else 0,
                    "RMSE": sklearn.metrics.mean_squared_error(pred_df.query(f"new_{method}_at_risk")["measurement"], pred_df.query(f"new_{method}_at_risk")[f"{method}_model_pred" if model_type=="full" else "restricted_model_pred"], squared=False) if pred_df[f"new_{method}_at_risk"].sum()>0 else 0,
                    "r2": sklearn.metrics.r2_score(pred_df.query(f"new_{method}_at_risk")["measurement"], pred_df.query(f"new_{method}_at_risk")[f"{method}_model_pred" if model_type=="full" else "restricted_model_pred"]) if pred_df[f"new_{method}_at_risk"].sum()>0 else 0,
                    "individuals_count": len(pred_df.query(f"new_{method}_at_risk==True and at_risk=={at_risk}"))
                })
res_df = pd.DataFrame.from_dict(res)

# %%
plot1 = (
    pn.ggplot(res_df.query("model_type=='full'"), pn.aes(y="individuals", x="method"))
    + pn.ggtitle(f"Number of individuals where prediction differs by more than {distance_std} standard deviation(s) from common PRS and measuremt is extreme (1%)")
    + pn.geom_bar(stat="identity")
    + pn.facet_wrap("phenotype", nrow=1)
    + pn.theme(figure_size=(16, 4))
)
display(plot1)

# %%
plot2 = (
    pn.ggplot(res_df, pn.aes(y="RMSE", x="method", fill="model_type"))
    + pn.ggtitle(f"RMSE of predictions where prediction differs by more than {distance_std} standard deviation(s) from common PRS")
    + pn.geom_bar(stat="identity", position="dodge")
    + pn.facet_wrap("phenotype", nrow=1)
    + pn.theme(figure_size=(12, 4))
)
display(plot2)

# %%
