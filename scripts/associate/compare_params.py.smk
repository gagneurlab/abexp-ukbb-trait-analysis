SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


OUTPUT_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}/fset={{feature_set}}'''

rule associate__compare_params:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        touch_file=touch(f"{OUTPUT_BASEPATH}/compare_params.done"),
    input:
        protein_coding_genes_pq=config["protein_coding_genes_pq"],
        associations_pq=f"{OUTPUT_BASEPATH}/associations.parquet",
        featureset_config=f"{OUTPUT_BASEPATH}/config.yaml",
        regression_done=f"{OUTPUT_BASEPATH}/done",
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


del OUTPUT_BASEPATH


