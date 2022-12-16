SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''

rule train_test_split:
    threads: 48
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        train_test_split_pq=f"{COVARIATES_BASEPATH}/train_test_split.parquet",
        train_test_split_tsv=f"{COVARIATES_BASEPATH}/train_test_split.tsv",
        train_test_split_pq_done=touch(f"{COVARIATES_BASEPATH}/train_test_split.parquet.done"),
    input:
        covariates_pq=f'''{COVARIATES_BASEPATH}/covariates.parquet''',
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


del COVARIATES_BASEPATH
