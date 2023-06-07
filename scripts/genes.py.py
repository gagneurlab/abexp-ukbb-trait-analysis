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
# # code autoreload
# # %load_ext autoreload
# # %autoreload 2
import os
import sys
import shutil

from pprint import pprint

import collections
import random
import math
import numpy as np
import pandas as pd

import json
import yaml

# import joblib

try:
    from tqdm.notebook import tqdm
except:
    from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib

import pyranges as pr

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

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
        rule_name = 'genes',
        default_wildcards={
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
os.getcwd()

# %% [markdown]
# # Load input data

# %% [markdown]
# ## Load protein coding genes

# %%
snakemake.input

# %%
gtf_df = (
    pr.read_gtf(snakemake.input["gtf_file"], as_df=True)
    .rename(columns={
        "biotype": "transcript_biotype",
        "transcript_type": "transcript_biotype",
    })
)
gtf_df

# %%
genes = gtf_df.query("Feature == 'gene'")
genes = genes.assign(gene_id=genes["gene_id"].str.split(".").apply(lambda s: s[0]))
genes

# %%
protein_coding_genes = genes.query("gene_type == 'protein_coding'")[["gene_id", "gene_name"]].drop_duplicates()
protein_coding_genes

# %%
# genes.query("gene_type == 'protein_coding'").groupby(["gene_id"]).size().rename("size").to_frame().query("size > 1").join(genes.set_index("gene_id"), how="left").set_index("Chromosome", append=True)


# %%
# protein_coding_transcripts = transcripts.query("transcript_biotype == 'protein_coding'")
# protein_coding_transcripts = protein_coding_transcripts.assign(gene_id=protein_coding_transcripts["gene_id"].str.split(".").apply(lambda s: s[0]))
# protein_coding_transcripts = protein_coding_transcripts.set_index("gene_id").sort_index()
# protein_coding_transcripts

# %%
# protein_coding_genes = protein_coding_transcripts.index.unique()
# protein_coding_genes

# %% [markdown]
# # Save output

# %%
snakemake.output

# %%
protein_coding_genes.to_csv(snakemake.output["protein_coding_genes_csv"], index=False)

# %%
protein_coding_genes.to_parquet(snakemake.output["protein_coding_genes_pq"], index=False)

# %%
