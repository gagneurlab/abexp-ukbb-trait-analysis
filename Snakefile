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
configfile: "config_v7.yaml"

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
UKBB_DECODED_PHENOTYPES_DIR = config["ukbb_decoded_phenotypes_dir"]

UKBB_SKIP_METADATA_UPDATE = bool(os.environ.get("UKBB_SKIP_METADATA_UPDATE", "false"))

phenotype_dirs, ukbb_codes = glob_wildcards(UKBB_RAW_PHENOTYPES_DIR + "/{phenotype_dir}/{ukbb_code}.tab")

        
# include: config['wBuildPath'] + "/wBuild.snakefile"
# include: 'additional_rules/rules/data_preprocessing.smk'

# associations = pd.read_csv('/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/signif_genes.tsv', sep = '\t')


include: 'scripts/__init__.smk'

hdl_cholesterol=expand(
    rules.associate__compare_genebass.output,
    phenotype_col="HDL_cholesterol", 
    feature_set=["LOFTEE_pLoF", "AbExp_all_tissues", "max_AbExp", "median_AbExp"],
    covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS", "sex_age_genPC_BMI_smoking_CLMP_PRS"],
)

hdl_cholesterol_term_pvals=expand(
    rules.associate__compare_params.output,
    phenotype_col="HDL_cholesterol", 
    feature_set=["AbExp_all_tissues",],
    covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS", "sex_age_genPC_BMI_smoking_CLMP_PRS"],
)

all_phenotypes = [
    # "Asthma",
    "BodyMassIndex",
    # "CAD_HARD",
    # "CAD_SOFT",
    # "c_reactive_protein", # missing data
    # "Diabetes",
    "glycated_haemoglobin_hba1c", # blood sugar -> diabetes risk
    "HDL_cholesterol",
    "LDL_direct",
    "Lipoprotein_A",
    # # "severe_LDL",
    "standing_height",
    "systolic_blood_pressure", # bad results
    "Triglycerides",
]
all_phenotypes_output = [
    *expand(
        rules.associate__compare_genebass.output,
        phenotype_col=all_phenotypes,
        feature_set=[
            "LOFTEE_pLoF",
            "AbExp_all_tissues",
            "max_AbExp",
            "median_AbExp"
        ],
        covariates=[
            "sex_age_genPC",
            "sex_age_genPC_CLMP",
            "sex_age_genPC_CLMP_PRS",
            "sex_age_genPC_BMI_smoking_CLMP_PRS",
            # "randomized_sex_age_genPC_CLMP_PRS",
        ],
    ),
    *expand(
        rules.associate__compare_params.output,
        phenotype_col=all_phenotypes,
        feature_set=["AbExp_all_tissues",],
        covariates=[
            "sex_age_genPC",
            "sex_age_genPC_CLMP",
            "sex_age_genPC_CLMP_PRS",
            "sex_age_genPC_BMI_smoking_CLMP_PRS",
            # "randomized_sex_age_genPC_CLMP_PRS",
        ],
    ),
    *expand(
        rules.associate__qq_plot.output,
        phenotype_col=all_phenotypes,
        feature_set=[
            "LOFTEE_pLoF",
            "AbExp_all_tissues",
            "max_AbExp",
            "median_AbExp"
        ],
        covariates=["randomized_sex_age_genPC_BMI_smoking_CLMP_PRS", ],
    ),
    # *expand(
    #     rules.associate__polygenic_risk_score.output,
    #     phenotype_col=all_phenotypes,
    #     feature_set=["LOFTEE_pLoF", "AbExp_all_tissues", "max_AbExp", "median_AbExp"],
    #     covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
    # )
]

rule all:
    input:
        expand(rules.read_phenotypes.output, zip, pheno_dir=phenotype_dirs, ukbb_code=ukbb_codes),
        # rules.merge_phenotype_metadata.output,
        expand(rules.filter_genebass.output, genebass_version=["300k", "500k"]),
#         rules.Index.output, 
        # *hdl_cholesterol,
        # *hdl_cholesterol_term_pvals,
        #*all_phenotypes_output,
        expand(rules.compare_associations.output, comparison=[
            # "all",
            "paper_figure",
            # "paper_figure_all_traits",
            "bmi_smoking_best_tissue_test",
        ]),
        expand(rules.compare_risk_scores.output, comparison=[
            # "all",
            "paper_figure",
        ]),
        # expand(
        #     rules.associate__polygenic_risk_score.output, 
        #     feature_set=[
        #         "LOFTEE_pLoF",
        #         "AbExp_all_tissues",
        #         "max_AbExp",
        #         "median_AbExp"
        #     ],
        #     phenotype_col=all_phenotypes,
        #     covariates=["sex_age_genPC_CLMP_PRS"]
        # ),
        rules.compose_paper_figure.output, 

localrules: all
