import os
import sys
import numpy as np
import pandas as pd
import glob
import yaml
import json
import pathlib

include: 'snakefile_utils.smk'

# import wbuild
# config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

workdir: "./"
configfile: "config.yaml"

# sync config directory
eprint("Syncing config directory...")
os.system(
    f'''rsync -r --links --partial --update --checksum "{config["config_dir"]}/" "{config["output_basedir"]}/"'''
)
eprint("Syncing config directory done!")

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
CONDA_ENV_YAML_DIR = f"{SNAKEMAKE_DIR}/envs"

# CACHE_DIR = config["dirs"]["CACHE_DIR"]
# RAW_DATA_DIR = config["dirs"]["RAW_DATA_DIR"]
# PROCESSED_DATA_DIR = config["dirs"]["PROCESSED_DATA_DIR"] # "/s/project/rep/processed"

UKBB_RAW_PHENOTYPES_DIR = config["ukbb_raw_phenotypes_dir"]
UKBB_PROCESSED_PHENOTYPES_DIR = config["ukbb_processed_phenotypes_dir"]




phenotype_dirs, ukbb_codes = glob_wildcards(UKBB_RAW_PHENOTYPES_DIR + "/{phenotype_dir}/{ukbb_code}.tab")

        
# include: config['wBuildPath'] + "/wBuild.snakefile"
# include: 'additional_rules/rules/data_preprocessing.smk'

# associations = pd.read_csv('/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/signif_genes.tsv', sep = '\t')


include: 'scripts/__init__.smk'

hdl_cholesterol=expand(
    rules.associate__compare_genebass.output,
    phenotype_col="hdl_cholesterol_f30760_0_0", 
    feature_set=["LOFTEE_pLoF", "AbExp_pivot", "max_AbExp", "median_AbExp"],
    covariates=["sex+age+genPC", "sex+age+genPC+CLMP"],
)

rule all:
    input:
        expand(rules.read_phenotypes.output, zip, pheno_dir=phenotype_dirs, ukbb_code=ukbb_codes),
        rules.merge_phenotype_metadata.output,
        rules.filter_genebass.output,
#         rules.Index.output, 
#         expand('/s/project/bayesRare/UKBB_Splicing_Analysis/results/{phenocode}/boxplot.png', phenocode = set(associations.phenocode))
        *hdl_cholesterol,


localrules: all