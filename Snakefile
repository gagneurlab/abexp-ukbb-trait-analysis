import os
import sys
import numpy as np
import pandas as pd
import glob
import yaml
import json
import pathlib

include: 'snakefile_utils.smk'

workdir: "./"
configfile: "config_v8.yaml"

# sync config directory
eprint("Syncing config directory...")
os.system(
    f'''rsync -r --links --partial --update --checksum "{config["config_dir"]}/" "{config["output_basedir"]}/"'''
)
eprint("Syncing config directory done!")

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
CONDA_ENV_YAML_DIR = f"{SNAKEMAKE_DIR}/envs"

UKBB_RAW_PHENOTYPES_DIR = config["ukbb_raw_phenotypes_dir"]
UKBB_PROCESSED_PHENOTYPES_DIR = config["ukbb_processed_phenotypes_dir"]
UKBB_DECODED_PHENOTYPES_DIR = config["ukbb_decoded_phenotypes_dir"]

UKBB_SKIP_METADATA_UPDATE = bool(os.environ.get("UKBB_SKIP_METADATA_UPDATE", "false"))

phenotype_dirs, ukbb_codes = glob_wildcards(UKBB_RAW_PHENOTYPES_DIR + "/{phenotype_dir}/{ukbb_code}.tab")


include: 'scripts/__init__.smk'

hdl_cholesterol=expand(
    rules.associate__compare_genebass.output,
    phenotype_col="HDL_cholesterol", 
    feature_set=["LOFTEE_pLoF", "AbExp_all_tissues", "max_AbExp", "median_AbExp"],
    covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
)

hdl_cholesterol_term_pvals=expand(
    rules.associate__compare_params.output,
    phenotype_col="HDL_cholesterol", 
    feature_set=["AbExp_all_tissues",],
    covariates=["sex_age_genPC", "sex_age_genPC_CLMP", "sex_age_genPC_CLMP_PRS"],
)

all_traits = [
    "Alanine_aminotransferase",
    "Albumin",
    "Alkaline_phosphatase",
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Aspartate_aminotransferase",
    "Calcium",
    "Cholesterol",
    "Creatinine",
    "Cystatin_C",
    "Direct_bilirubin",
    "Eosinophill_count",
    "Erythrocyte_distribution_width",
    "Gamma_glutamyltransferase",
    "Glucose",
    "HDL_cholesterol",
    "Haematocrit_percentage",
    "IGF1",
    "LDL_direct",
    "Leukocyte_count",
    "Lipoprotein_A",
    "Lymphocyte_percentage",
    "Mean_corpuscular_haemoglobin",
    "Mean_corpuscular_volume",
    "Mean_reticulocyte_volume",
    "Mean_sphered_cell_volume",
    "Monocyte_count",
    "Neutrophill_percentage",
    "Phosphate",
    "Platelet_count",
    "Reticulocyte_count",
    "SHBG",
    "Testosterone",
    "Thrombocyte_volume",
    "Total_bilirubin",
    "Triglycerides",
    "Urate",
    "Vitamin_D",
    "c_reactive_protein",
    "glycated_haemoglobin_hba1c",
]

rule all:
    input:
        expand(rules.read_phenotypes.output, zip, pheno_dir=phenotype_dirs, ukbb_code=ukbb_codes),
        # rules.merge_phenotype_metadata.output,
        expand(rules.filter_genebass.output, genebass_version=["300k", "500k"]),
        # *hdl_cholesterol,
        # *hdl_cholesterol_term_pvals,
        expand(rules.compare_associations.output, comparison=[
            "paper_figure",
            "paper_figure_randomized",
        ]),
        expand(rules.compare_risk_scores.output, comparison=[
            "paper_figure",
        ]),
        expand(
            rules.associate__qq_plot.output,
            phenotype_col=all_traits,
            feature_set=[
                "LOFTEE_pLoF",
                "AbExp_all_tissues",
                "max_AbExp",
                "median_AbExp"
            ],
            covariates=["randomized_sex_age_genPC_CLMP_PRS", ],
        ),
        expand(
            rules.compare_stringdb_image.output,
            phenotype=all_traits,
            comparison="AbExp_all_tissues_vs_LOFTEE",
            covariates=["sex_age_genPC_CLMP_PRS", ],
        ),
        # expand(
        #     rules.associate__polygenic_risk_score.output, 
        #     feature_set=[
        #         "LOFTEE_pLoF",
        #         "AbExp_all_tissues",
        #         "max_AbExp",
        #         "median_AbExp"
        #     ],
        #     phenotype_col=all_traits,
        #     covariates=["sex_age_genPC_CLMP_PRS"]
        # ),
        rules.compose_paper_figure.output, 

localrules: all
