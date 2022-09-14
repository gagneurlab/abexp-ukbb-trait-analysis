SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


OUTPUT_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''
TEMPLATE_FILE=f"{SNAKEFILE_DIR}/covariates@{{covariates}}.yaml"


rule covariates:
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        metadata_pq=f"{OUTPUT_BASEPATH}/metadata.parquet",
        metadata_tsv=f"{OUTPUT_BASEPATH}/metadata.tsv",
        metadata_pq_done=touch(f"{OUTPUT_BASEPATH}/metadata.parquet.done"),
        covariates_config=f"{OUTPUT_BASEPATH}/config.yaml",
        covariates_pq=f"{OUTPUT_BASEPATH}/covariates.parquet",
        clumping_variants_pq=f"{OUTPUT_BASEPATH}/clumping_variants.parquet",
    input:
        samples_txt=config["samples"],
        phenotype_metadata_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
        covariates_config=TEMPLATE_FILE,
        # clumping
        mac_index_variants=config["mac_index_variants"],
        genome_annotation=config["gtf_file"],
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


del OUTPUT_BASEPATH
del TEMPLATE_FILE
