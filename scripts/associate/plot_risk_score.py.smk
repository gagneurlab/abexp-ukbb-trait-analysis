SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''
ASSOCIATION_BASEPATH=f'''{COVARIATES_BASEPATH}/fset={{feature_set}}'''

OUTPUT_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}/fset={{feature_set}}/polygenic_risk_score'''


rule associate__plot_risk_score:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    input: **rules.associate__polygenic_risk_score.output
    output:
        plotting_done=touch(f"{OUTPUT_BASEPATH}/plotting.done"),
        # protein_coding_genes_pq=config["protein_coding_genes_pq"],
        # associations_pq=f"{ASSOCIATION_BASEPATH}/associations.parquet",
        # featureset_config=f"{ASSOCIATION_BASEPATH}/config.yaml",
        # regression_done=f"{ASSOCIATION_BASEPATH}/done",
        # # covariates
        # covariates_pq=f'''{COVARIATES_BASEPATH}/covariates.parquet''',
        # # clumping
        # clumping_variants_pq=f'''{COVARIATES_BASEPATH}/clumping_variants.parquet''',
        # # sample splits
        # samples_pq=f"{COVARIATES_BASEPATH}/samples.parquet",
        # samples_pq_done=f"{COVARIATES_BASEPATH}/samples.parquet.done",
        # train_test_split_pq=f"{COVARIATES_BASEPATH}/train_test_split.parquet",
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        output_basedir=OUTPUT_BASEPATH,
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


del COVARIATES_BASEPATH
del ASSOCIATION_BASEPATH
del OUTPUT_BASEPATH

