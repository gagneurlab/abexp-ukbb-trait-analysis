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

import polars as pl


# %%
def str_presenter(dumper, data):
    """configures yaml for dumping multiline strings
    Ref: https://stackoverflow.com/questions/8640959/how-can-i-control-what-scalar-form-pyyaml-uses-for-my-data"""
    if data.count('\n') > 0:  # check for multiline string
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='|')
    return dumper.represent_scalar('tag:yaml.org,2002:str', data)

yaml.add_representer(str, str_presenter)

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
        rule_name = 'covariates',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            "covariates": "sex+age+genPC+CLMP",
        }
    )

# %%
print(json.dumps(snakemake.__dict__, indent=2, default=str))


# %% [markdown]
# # Load configuration

# %%
def recursive_format(data, params, fail_on_unknown=False):
    try:
        if isinstance(data, str):
            return data.format(**params)
        elif isinstance(data, dict):
            return {k: recursive_format(v, params) for k, v in data.items()}
        elif isinstance(data, list):
            return [recursive_format(v, params) for v in data]
        else:
            if fail_on_unknown:
                raise ValueError("Handling of data type not implemented: %s" % type(data))
            else:
                return data
    except Exception as e:
        print(data)
        print(params)
        raise e

class AttrDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__



# %%
with open(snakemake.input["covariates_config"], "r") as fd:
    config = yaml.safe_load(fd)

config = recursive_format(config, params=dict(params=AttrDict(snakemake.params), wildcards=AttrDict(snakemake.wildcards), config=AttrDict(snakemake.config)))

# %%
print(json.dumps(config, indent=2, default=str))

# %% [markdown]
# ## set defaults

# %%
config["add_clumping"] = config.get("add_clumping", False)
config["add_clumping"]

# %%
config["clumping_gene_padding"] = config.get("clumping_gene_padding", 0)
config["clumping_gene_padding"]

# %%
config["add_PRS"] = config.get("add_PRS", False)
config["add_PRS"]

# %%
restricted_formula_template = config["restricted_formula"]
restricted_formula = "\n + ".join(restricted_formula_template)

config["restricted_formula"] = restricted_formula
print(restricted_formula)

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %%
import re

m = re.search('.*_f(.+?)_.*', phenotype_col)
if m:
    phenocode = m.group(1)
else:
    raise ValueError("Cannot find phenocode!")
phenocode

# %%
config["phenocode"] = phenocode

# %%
print(json.dumps(config, indent=2, default=str))

# %%
with open(snakemake.output["covariates_config"], "w") as fd:
    yaml.dump(config, fd, sort_keys=False)

# %% [markdown] {"jp-MarkdownHeadingCollapsed": true, "tags": []}
# # Read metadata

# %%
phenotype_metadata_df = pd.read_parquet(snakemake.input["phenotype_metadata_pq"])
phenotype_metadata_df = phenotype_metadata_df.rename(columns={"field.showcase": "phenocode"})
phenotype_metadata_df

# %%
import patsy

model_desc = patsy.ModelDesc.from_formula(restricted_formula)

covariate_cols = list(dict.fromkeys(
    [factor.name() for term in model_desc.rhs_termlist for factor in term.factors]
))
covariate_cols

# %%
meta_cols = [c for c in covariate_cols if c not in {
    "{CLMP}",
    "{PRS}",
}]
phenotype_metadata_subset = phenotype_metadata_df.set_index("col.name").loc[[
    phenotype_col,
    *meta_cols,
#     *phenotype_column,
]]
phenotype_metadata_subset

# %%
phenotype_metadata_subset.reset_index().to_parquet(snakemake.output["metadata_pq"], index=False)

# %%
phenotype_metadata_subset.reset_index().to_csv(snakemake.output["metadata_tsv"], index=False, header=True, sep="\t")

# %% [markdown]
# # Generate covariates

# %% [markdown]
# ## read samples

# %%
with open(snakemake.input["samples_txt"], "r") as fd:
    samples = fd.read().splitlines()
samples[:10]

# %%
samples_df = pl.DataFrame({
    "individual": sorted(samples)
})
samples_df

# %% [markdown] {"tags": []}
# ## phenotypes

# %%
data_dfs = []
for data_path, group_df in phenotype_metadata_subset.groupby("data_path"):
    data_df = (
        pl.scan_parquet(data_path)
        .select([
            "eid",
            *group_df.index.to_list()
        ])
    )
    data_dfs.append(data_df)

# %%
len(data_dfs)

# %%
import functools
data_df = functools.reduce(lambda df1, df2: df1.join(df2, on="eid", how="outer"), data_dfs)

# %%
# remove NA values from the phenotype column
data_df = data_df.filter(pl.col(phenotype_col).is_not_null())

# %%
# rename eid col
data_df = data_df.rename({"eid": "individual"})

# %%
# make sure that sample_id is a string
data_df = data_df.with_column(pl.col("individual").cast(pl.Utf8).alias("individual"))

# %%
# change order of columns
data_df = data_df.select([
    "individual",
    *phenotype_metadata_subset.index
])

# %%
data_df.schema

# %% [markdown] {"tags": []}
# ## clumping

# %% [markdown]
# ### read clumping variant calls

# %% {"tags": []}
mac_index_vars = (
    pl.scan_parquet(snakemake.input["mac_index_variants"])
    .rename({"ID": "individual"})
)
mac_index_vars.schema

# %% [markdown]
# ### parse clumping variants

# %%
vcf_pattern = r"(.+):(\d+):(.*)>(.*)"
variants = (
    pl.DataFrame({"variant": pl.Series(mac_index_vars.columns)})
    .filter(pl.col("variant") != pl.lit("individual"))
    .with_columns([
        pl.col("variant").str.extract(vcf_pattern, 1).alias("Chromosome"),
        pl.col("variant").str.extract(vcf_pattern, 2).alias("pos").cast(pl.Int32),
        pl.col("variant").str.extract(vcf_pattern, 3).alias("ref"),
        pl.col("variant").str.extract(vcf_pattern, 4).alias("alt"),
    ])
    .with_columns([
        (pl.col("pos") - 1).alias("Start"),
        (pl.col("pos") - 1 + pl.col("ref").str.lengths()).alias("End"),
    ])
)
variants

# %%
import pyranges as pr

# %%
genome_annotation = pr.read_gtf(snakemake.input["genome_annotation"], as_df=True)
gene_annoation = genome_annotation[genome_annotation.Feature == "gene"]
gene_annoation

# %%
padding = config["clumping_gene_padding"]
padding

# %%
gene_annotation_expanded = gene_annoation[["Chromosome", "Start", "End", "gene_id"]].copy()
gene_annotation_expanded["gene"] = gene_annotation_expanded["gene_id"].str.split(".").apply(lambda s: s[0])
gene_annotation_expanded["gene_start"] = gene_annotation_expanded.Start
gene_annotation_expanded["gene_end"] = gene_annotation_expanded.End
gene_annotation_expanded.Start = np.fmax(0, gene_annotation_expanded.Start - padding)
gene_annotation_expanded.End = gene_annotation_expanded.End + padding
gene_annotation_expanded

# %%
clumping_meta_df = (
    pr.PyRanges(variants.to_pandas())
    .join(pr.PyRanges(gene_annotation_expanded), suffix="_gene")
    .drop(["End_gene", "Start_gene"])
    .df
    .loc[:, ["variant", "Chromosome", "gene", "gene_start", "gene_end"]]
)
clumping_meta_df

# %% [markdown]
# ### write clumping variants

# %%
clumping_meta_df.to_parquet(snakemake.output["clumping_variants_pq"])

# %% [markdown]
# ## join everything

# %% {"tags": []}
covariates_df = (
    samples_df.lazy()
    .join(data_df, on="individual", how="left")
    .join(mac_index_vars, on="individual", how="left")
    .filter(pl.col(phenotype_col).is_not_null())
    .collect()
)
covariates_df

# %%
covariates_df.write_parquet(
    snakemake.output["covariates_pq"],
    compression="snappy",
    use_pyarrow=True,
)

# %%
