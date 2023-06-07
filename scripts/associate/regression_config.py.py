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
            "covariates": "sex+age+genPC+CLMP",
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
## add covariate config
with open(snakemake.input["covariate_config"], "r") as fd:
    covariate_config = yaml.safe_load(fd)

# %%
config["covariates"] = covariate_config

# %%
print(json.dumps(config, indent=2, default=str))

# %%
with open(snakemake.output["config"], "w") as fd:
    yaml.dump(config, fd, sort_keys=False)

# %%
