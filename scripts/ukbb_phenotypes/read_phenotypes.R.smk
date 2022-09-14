SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

rule read_phenotypes:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        output_data_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/{{pheno_dir}}/{{ukbb_code}}.data.parquet",
        output_meta_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/{{pheno_dir}}/{{ukbb_code}}.meta.parquet",
    input:
        input_data_dir=f"{UKBB_RAW_PHENOTYPES_DIR}/{{pheno_dir}}/",
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
    wildcard_constraints:
        ds_dir="[^/]+",
    conda: f'{CONDA_ENV_YAML_DIR}/ukbb-trait-analysis-R.yaml'
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    shell:
        "Rscript {params.nb_script}.R {input.input_data_dir} {wildcards.ukbb_code} {output.output_data_pq} {output.output_meta_pq}"