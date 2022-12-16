SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''

rule samples:
    threads: 48
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        samples_pq=f"{COVARIATES_BASEPATH}/samples.parquet",
        samples_tsv=f"{COVARIATES_BASEPATH}/samples.tsv",
        samples_pq_done=touch(f"{COVARIATES_BASEPATH}/samples.parquet.done"),
    input:
        samples_txt=config["samples"],
        phenotype_metadata_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
        # covariates_config=f"{COVARIATES_BASEPATH}/config.yaml",
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


del COVARIATES_BASEPATH
