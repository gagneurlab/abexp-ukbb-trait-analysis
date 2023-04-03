SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''
TEMPLATE_FILE=f"{SNAKEFILE_DIR}/covariates@{{covariates}}.yaml"


def _covariates_input_fn(wildcards, yaml_path=TEMPLATE_FILE, covariates_basepath=COVARIATES_BASEPATH):
    import yaml
    
    yaml_path = yaml_path.format(covariates=wildcards.covariates)
    
    try:
        with open(yaml_path, "r") as fd:
            covariates_config=yaml.safe_load(fd)
    except:
        covariates_config = {}
    
    input_files=dict(
        samples_pq=f"{covariates_basepath}/samples.parquet",
        samples_pq_done=f"{covariates_basepath}/samples.parquet.done",
        decoded_phenotype_pq=f"{UKBB_DECODED_PHENOTYPES_DIR}/{{phenotype_col}}/data.parquet",
        phenotype_metadata_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
        covariates_config=ancient(yaml_path),
        # clumping
        # mac_index_variants=config["mac_index_variants"],
        clumping_genome_annotation=config["clumping_gtf_file"],
        # PRS scores
        PRS_score_mapping_csv=ancient(config["PRS_score_mapping"]),
        register_phenotypes=[],
    )
    
    if "register_phenotypes" in covariates_config:
        input_files["register_phenotypes"] = [
            f"{UKBB_DECODED_PHENOTYPES_DIR}/{phenotype}/data.parquet" for phenotype in covariates_config["register_phenotypes"]
        ]
    
    return input_files


rule covariates:
    threads: 48
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        metadata_pq=f"{COVARIATES_BASEPATH}/metadata.parquet",
        metadata_tsv=f"{COVARIATES_BASEPATH}/metadata.tsv",
        metadata_pq_done=touch(f"{COVARIATES_BASEPATH}/metadata.parquet.done"),
        covariates_config=f"{COVARIATES_BASEPATH}/config.yaml",
        covariates_pq=f"{COVARIATES_BASEPATH}/covariates.parquet",
        clumping_variants_pq=f"{COVARIATES_BASEPATH}/clumping_variants.parquet",
    input: unpack(_covariates_input_fn)
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        # age_col='age_when_attended_assessment_centre_f21003_0_0',
        # sex_col='sex_f31_0_0',
        # year_of_birth_col='year_of_birth_f34_0_0',
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"

        
rule covariates_to_ipc:
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        covariates_ipc=temp(f"{COVARIATES_BASEPATH}/covariates.feather"),
    input:
        covariates_pq=f"{COVARIATES_BASEPATH}/covariates.parquet",
    wildcard_constraints:
        covariates="[^/]+",
    run:
        import polars as pl
        (
            pl.read_parquet(input["covariates_pq"])
            .write_ipc(output["covariates_ipc"])
        )

del COVARIATES_BASEPATH
del TEMPLATE_FILE
del _covariates_input_fn
