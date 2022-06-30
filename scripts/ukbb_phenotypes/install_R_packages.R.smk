SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

rule read_phenotypes__setup_R_env:
    threads: 2
    resources:
        ntasks=1,
        mem_mb=4000,
    output:
        install_R_packages_done=touch(INSTALL_R_PACKAGES_DONE),
    input:
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
    wildcard_constraints:
    conda: f'{CONDA_ENV_YAML_DIR}/ukbb-trait-analysis-R.yaml'
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    shell:
        "Rscript {params.nb_script}.R"

