SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


rule filter_genebass:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        results_filtered_pq=directory(config["genebass_results_filtered_pq"]),
        filtered_phenotypes_pq=config["genebass_filtered_phenotypes_pq"],
    input:
        results_pq=config["genebass_results_pq"],
        phenotype_metadata_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
    # conda: f'{CONDA_ENV_YAML_DIR}/ukbb-trait-analysis-py.yaml'
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        bonferroni_cutoff=float(config["genebass_bonferroni_cutoff"])
    wildcard_constraints:
        ds_dir="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"
        
        
