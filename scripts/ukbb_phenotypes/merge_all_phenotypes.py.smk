SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

rule merge_all_phenotypes:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        latest_data_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.data.parquet",
    input:
        latest_meta_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
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