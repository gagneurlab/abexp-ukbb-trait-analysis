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

# %%
import os
import numpy as np
import pandas as pd

import json
import yaml

import re

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
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'covariates',
        default_wildcards={
            # "phenotype_col": "Asthma",
            # "phenotype_col": "Lipoprotein_A",
            # "phenotype_col": "Diabetes",
            "phenotype_col": "HDL_cholesterol",
            # "phenotype_col": "Triglycerides",
            # "phenotype_col": "standing_height",
            # "phenotype_col": "systolic_blood_pressure_f4080_0_0",
            # "covariates": "randomized_sex_age_genPC_CLMP_PRS",
            # "covariates": "randomizedOLS_sex_age_genPC_CLMP_PRS",
            # "covariates": "sex_age_genPC_CLMP_PRS",
            "covariates": "sex_age_genPC_BMI_smoking_CLMP_PRS",
            # "covariates": "sex+age+genPC+CLMP",
            # "covariates": "sex_age_genPC",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))


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
config["clumping_use_pca"] = config.get("clumping_use_pca", False)
config["clumping_use_pca"]

# %%
config["add_PRS"] = config.get("add_PRS", False)
config["add_PRS"]

# %%
config["randomize_phenotype"] = config.get("randomize_phenotype", False)
config["randomize_phenotype"]

# %%
phenotype_col = snakemake.wildcards["phenotype_col"]
phenotype_col

# %%
config["force_quantitative_regression"] = config.get("force_quantitative_regression", False)
config["force_quantitative_regression"]

# %% [markdown]
# ## PRS / Clumping / phenocode

# %%
prs_score_mapping = (
    pl.read_csv(snakemake.input["PRS_score_mapping_csv"])
    .filter(pl.col("phenotype") == pl.lit(phenotype_col))
)
prs_score_mapping

# %%
mac_index_variants_path = prs_score_mapping.select("clumping_var_file_path").unique().to_numpy().item()
mac_index_variants_path

# %%
phenocode = str(prs_score_mapping.select("genebass_phenocode").unique().to_numpy().item())
phenocode

# %%
config["phenocode"] = phenocode

# %% [markdown]
# ## restricted formula

# %%
restricted_formula_template = [
    *config["restricted_formula"],
    *(prs_score_mapping["pgs_id"] if config["add_PRS"] else []),
]
restricted_formula = "\n + ".join(restricted_formula_template)

config["restricted_formula"] = restricted_formula
print(restricted_formula)

# %%
print(json.dumps(config, indent=2, default=str))

# %% [markdown]
# ## write config.yaml

# %%
with open(snakemake.output["covariates_config"], "w") as fd:
    yaml.dump(config, fd, sort_keys=False)

# %% [markdown]
# # Read metadata

# %% [markdown]
# ## phenotypes

# %%
phenotype_metadata_df = pd.read_parquet(snakemake.input["phenotype_metadata_pq"])
phenotype_metadata_df = phenotype_metadata_df.rename(columns={"field.showcase": "phenocode"})
phenotype_metadata_df

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
covariate_cols = get_variables_from_formula(restricted_formula, lhs=False)
covariate_cols

# %%
meta_cols = [
    c for c in covariate_cols 
    if c not in prs_score_mapping["pgs_id"].to_list() 
    and c not in config.get("register_phenotypes", [])
]
phenotype_metadata_subset = phenotype_metadata_df.set_index("col.name").loc[[
    # phenotype_col,
    *meta_cols,
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
samples_df = pl.read_parquet(snakemake.input["samples_pq"]).sort("individual")
samples_df

# %% [markdown]
# ## phenotypes

# %%
phenotype_df = pl.read_parquet(snakemake.input["decoded_phenotype_pq"]).sort("eid").lazy()
phenotype_df.schema

# %%
data_dfs = [phenotype_df]

for data_path in snakemake.input["register_phenotypes"]:
    data_df = pl.scan_parquet(data_path)
    print(data_df.schema)
    data_dfs.append(data_df)

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
data_df = data_df.with_columns(pl.col("individual").cast(pl.Utf8).alias("individual"))

# %%
# shuffle phenotype column if requested
if config["randomize_phenotype"]:
    print("shuffling phenotype...")
    data_df = data_df.with_columns(pl.col(phenotype_col).shuffle(seed=42).alias(phenotype_col))

# %%
# shuffle phenotype column if requested
if config["force_quantitative_regression"]:
    print("forcing quantitative phenotype...")
    data_df = data_df.with_columns(pl.col(phenotype_col).cast(pl.Float32()).alias(phenotype_col))

# %%
# change order of columns
data_df = data_df.select([
    "individual",
    snakemake.wildcards["phenotype_col"],
    *[c for c in covariate_cols if c in data_df]
])

# %%
data_df.schema

# %%
data_df = data_df.collect()

# %% [raw]
# ### quantile normalization

# %% [raw]
# import scipy.stats
#
# def irnt(col: pl.Series) -> np.array:
#     """
#     inverse rank normal transformation as used by PHESANT:
#     https://github.com/MRCIEU/PHESANT/blob/3f4a65d7fe93aaf01f3a4a3f39843562612a8d65/WAS/testContinuous.r#L243-L249
#     
#     :param x: Pandas Series with the (continuous) input values
#     :returns: irnt values
#     """
#     if col.is_float():
#         numPhenos = (col.is_not_nan() & col.is_not_null()).sum()
#     else:
#         numPhenos = col.is_not_null().sum()
#     
#     quantilePheno = (col.rank() - 0.5) / numPhenos
#     phenoIRNT = pl.Series(scipy.stats.norm.ppf(quantilePheno.to_numpy()))
#     return phenoIRNT

# %% [raw]
# columns_to_normalize = config.get("quantile_norm", [])
# columns_to_normalize

# %% [raw]
# data_df = data_df.with_columns([
#     irnt(data_df[col]).alias(col) for col in columns_to_normalize
# ])

# %%

# %%
data_df = data_df.lazy()

# %% [markdown]
# ## clumping

# %% [markdown]
# ### read clumping variant calls

# %%
mac_index_variants_path

# %%
mac_index_vars = (
    pl.scan_parquet(mac_index_variants_path)
    .rename({"IID": "individual"})
    .with_columns(pl.col("individual").str.extract("([^_]+).*", 1))
)
mac_index_vars.schema

# %%
mac_index_vars.select(mac_index_vars.columns[:6]).head().collect()

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
    .filter((
        pl.col("Chromosome").is_not_null()
        & pl.col("pos").is_not_null()
        & pl.col("ref").is_not_null()
        & pl.col("alt").is_not_null()
    ))
)
variants

# %%
# cast data to save memory
mac_index_vars = (
    mac_index_vars
    .select([
        "individual",
        *[(
            pl.col(c)
            .cast(pl.Int8)
            .fill_null(pl.lit(0))
            .alias(c)
        ) for c in variants["variant"]],
    ])
).collect().lazy()

# %%
import pyranges as pr

# %%
genome_annotation = pr.read_gtf(snakemake.input["clumping_genome_annotation"], as_df=True)
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
# ## PRS

# %%
pgs_dfs = []
for pgs_id, pgs_score_file_path in prs_score_mapping.select(["pgs_id", "pgs_score_file_path"]).rows():
    print(f"Loading '{pgs_id}'...")
    pgs_df = (
        pl.read_csv(pgs_score_file_path, separator="\t")
        .lazy()
        .rename({
            "#IID": "individual",
            "SCORE1_AVG": pgs_id,
        })
        .with_columns(pl.col("individual").str.extract("([^_]+).*", 1))
        .select([
            "individual",
            pgs_id,
        ])
    )
    pgs_dfs.append(pgs_df)

# %%
len(pgs_dfs)

# %%
import functools
pgs_df = functools.reduce(lambda df1, df2: df1.join(df2, on="individual", how="outer"), pgs_dfs)

# %% [markdown]
# # join everything

# %%
covariates_df = (
    samples_df.lazy()
    .join(data_df, on="individual", how="left")
    .join(pgs_df, on="individual", how="left")
    .join(mac_index_vars, on="individual", how="left")
    .filter(pl.col(phenotype_col).is_not_null())
    .fill_null(0)
    .collect()
)
covariates_df.schema

# %%
snakemake.output

# %%
(
    covariates_df
    .write_parquet(
        snakemake.output["covariates_pq"],
        compression="snappy",
        statistics=True,
        use_pyarrow=True,
    )
)

# %%
# (
#    pl.read_parquet(snakemake.output["covariates_pq"])
#    .write_ipc(
#        snakemake.output["covariates_pq"].replace("parquet", "feather"),
#    )
# )

# %%

