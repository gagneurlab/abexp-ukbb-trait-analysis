import os
import sys
import glob
import yaml
import json
import pathlib


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# from snakemk_util import recursive_format
def recursive_format(data, params, fail_on_unknown=False):
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


checkpoint file_depends:
    output: 
        file=touch("{file}.depends"),
    input:
        file="{file}",
    shell:
        "echo '{wildcards.file}'"


def require(file):
    """
    Makes sure that a certain input file exists before the rest of the rule is being executed.
    Technically, returns a function that takes wildcards as arguments.
    The resulting dependency adds ".depends" to the file path.
    Usage example:
        ```
        def read_config():
            [...]
        
        rule test:
            input:
                config_depends = require("config.yaml") # will effectively require the file "config.yaml.depends"
                config = lambda wildcards: read_config("config.yaml")
        ```
    This script does not fail, as `config_depends` has to exist before `config` is being executed.
    """
    return lambda wildcards: checkpoints.file_depends.get(file=file).output.file


def read_yaml_input(file, wildcards):
    """
    Reads the "input" section from a config yaml
    """
    import yaml
    with open(file, "r") as fd:
        config = yaml.safe_load(fd)

    retval = dict()
    if "snakemake" in config:
        if "input" in config["snakemake"]:
            retval = config["snakemake"]["input"]
            retval = recursive_format(retval, wildcards)
            
            # make sure that the output is a dictionary
            if isinstance(retval, list):
                retval = dict(yaml_input=retval)
    
    return retval


def require_yaml_input(file):
    """
    Reads the "input" section from a config yaml and adds it as additional input rules.
    Also, makes sure that the config file exists beforehand.
    See also `require`
    """
    def retval(wildcards):
        file_fmt = file.format(**wildcards)
        # make sure that yaml file exists:
        config_depends = checkpoints.file_depends.get(file=file_fmt).output

        return read_yaml_input(file_fmt, wildcards)
    return retval



localrules: file_depends
