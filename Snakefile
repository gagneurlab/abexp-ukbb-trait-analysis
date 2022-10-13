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
configfile: "config_v2.yaml"

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
    feature_set=["LOFTEE_pLoF", "AbExp_all_tissues", "max_AbExp", "median_AbExp"],
    covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
)

hdl_cholesterol_term_pvals=expand(
    rules.associate__compare_params.output,
    phenotype_col="hdl_cholesterol_f30760_0_0", 
    feature_set=["AbExp_all_tissues",],
    covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
)

all_phenotypes = [
    "hdl_cholesterol_f30760_0_0",
    "ldl_direct_f30780_0_0",
    # "c_reactive_protein_f30710_0_0",
    "triglycerides_f30870_0_0",
    "standing_height_f50_0_0",
    "body_mass_index_bmi_f21001_0_0",
    "systolic_blood_pressure_automated_reading_f4080_0_0",
]
all_phenotypes_output = [
    *expand(
        rules.associate__compare_genebass.output,
        phenotype_col=all_phenotypes,
        feature_set=["LOFTEE_pLoF", "AbExp_all_tissues", "max_AbExp", "median_AbExp"],
        covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
    ),
    *expand(
        rules.associate__compare_params.output,
        phenotype_col=all_phenotypes,
        feature_set=["AbExp_all_tissues",],
        covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
    ),
    *expand(
        rules.associate__polygenic_risk_score.output,
        phenotype_col=all_phenotypes,
        feature_set=["LOFTEE_pLoF", "AbExp_all_tissues", "max_AbExp", "median_AbExp"],
        covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
    )
]

rule all:
    input:
        expand(rules.read_phenotypes.output, zip, pheno_dir=phenotype_dirs, ukbb_code=ukbb_codes),
        rules.merge_phenotype_metadata.output,
        rules.filter_genebass.output,
#         rules.Index.output, 
#         expand('/s/project/bayesRare/UKBB_Splicing_Analysis/results/{phenocode}/boxplot.png', phenocode = set(associations.phenocode))
        *hdl_cholesterol,
        *hdl_cholesterol_term_pvals,
        *all_phenotypes_output,
        expand(rules.compare_associations.output, comparison=["all"]),


localrules: all
