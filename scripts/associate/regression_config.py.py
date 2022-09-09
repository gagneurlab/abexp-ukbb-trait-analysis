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
        rule_name = 'associate__regression_config',
        default_wildcards={
            "phenotype_col": "hdl_cholesterol_f30760_0_0",
            "feature_set": "LOFTEE_pLoF",
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
with open(snakemake.input["config_template"], "r") as fd:
    config = yaml.safe_load(fd)
config = recursive_format(config, params=dict(params=AttrDict(snakemake.params), wildcards=AttrDict(snakemake.wildcards), config=AttrDict(snakemake.config)))

# %%
config["snakemake"] = {
    "input": {
        "features": [f + "/data.parquet.done" for f in config["feature_sets"].values()]
    }
}

# %%
restricted_formula_template = []
for c in config["restricted_formula"]:
    if c == "{{clump}}":
        raise ValueError("Not yet implemented")
    elif c == "{{PRS}}":
        raise ValueError("Not yet implemented")
    else:
        restricted_formula_template.append(c)
restricted_formula = "\n + ".join(restricted_formula_template)

config["restricted_formula"] = restricted_formula
print(restricted_formula)

# %%
import patsy

model_desc = patsy.ModelDesc.from_formula(restricted_formula)

covariate_cols = list(dict.fromkeys(
    [factor.name() for term in model_desc.rhs_termlist for factor in term.factors]
))
covariate_cols

# %%
config["covariate_cols"] = covariate_cols

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
with open(snakemake.output["config"], "w") as fd:
    yaml.dump(config, fd)

# %%


# %%
