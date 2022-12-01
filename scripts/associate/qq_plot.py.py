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
        rule_name = 'associate__qq_plot',
        default_wildcards={
            "phenotype_col": "Triglycerides",
            # "phenotype_col": "triglycerides_f30870_0_0",
            # "phenotype_col": "standing_height_f50_0_0",
            # "phenotype_col": "body_mass_index_bmi_f21001_0_0",
            # "phenotype_col": "systolic_blood_pressure_automated_reading_f4080_0_0",
            # "phenotype_col": "hdl_cholesterol_f30760_0_0",
            # "feature_set": "LOFTEE_pLoF",
            "feature_set": "AbExp_all_tissues",
            # "feature_set": "LOFTEE_pLoF",
            "covariates": "randomized_sex_age_genPC_CLMP_PRS",
            # "covariates": "randomized_sex_age_genPC_CLMP",
            # "covariates": "randomized_sex_age_genPC",
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
regression_results_df = pl.scan_parquet(snakemake.input["associations_pq"] + "/*.parquet")

# %% [markdown]
# # Q-Q Plot

# %%
plot_df = regression_results_df.select(pl.col("lr_pval")).collect()

# %%
plot_df["lr_pval"].min()

# %%
snakemake.wildcards

# %%
plot = (
    pn.ggplot(plot_df, pn.aes(sample="lr_pval"))
    + pn.geom_abline(slope=1, linetype="dashed", color="red")
    + pn.stat_qq(distribution="uniform")
    + pn.scale_x_log10(limits=(10**-20, 1))
    + pn.scale_y_log10(limits=(10**-20, 1))
    + pn.labs(
        title="\n".join([
            f"""Q-Q plot of randomized p-values vs. random uniform distribution""",
            f"""phenotype: '{snakemake.wildcards["phenotype_col"]}'""",
            f"""feature set: '{snakemake.wildcards["feature_set"]}'""",
            f"""covariates: '{snakemake.wildcards["covariates"]}'""",
        ])
    )
    + pn.theme(title=pn.element_text(linespacing=1.4))
)
display(plot)

# %%
snakemake.output

# %%
pn.ggsave(plot, snakemake.output['qq_plot_pdf'], dpi=DPI)
pn.ggsave(plot, snakemake.output['qq_plot_png'], dpi=DPI)

# %% [raw]
# from IPython.display import Image
# display(Image(snakemake.params["output_basedir"] + "/phenotype_correlation.png"))

# %%
