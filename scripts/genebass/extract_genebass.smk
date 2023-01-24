SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

# RESULTS_MT = "/s/project/rep/processed/genebass/results.mt"
# RESULTS_TSV = "/s/project/rep/processed/genebass/results.tsv"
# RESULTS_PQ = "/s/project/rep/processed/genebass/results.parquet"

rule extract_genebass:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        results_pq=directory(config["genebass_results_pq"]),
    input:
        results_mt=config["genebass_results_mt"],
    conda: f'{CONDA_ENV_YAML_DIR}/ukbb-trait-analysis-hail.yaml'
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
    wildcard_constraints:
        ds_dir="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"
        

