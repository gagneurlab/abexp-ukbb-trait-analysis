SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


OUTPUT_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}/fset={{feature_set}}'''

rule associate__qq_plot:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        qq_plot_pdf=f"{OUTPUT_BASEPATH}/qq_plot.pdf",
        qq_plot_png=f"{OUTPUT_BASEPATH}/qq_plot.png",
    input:
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
