SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


rule associate:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        associations_pq=config["trait_associations"] + "/associations.parquet",
        data_pq=config["trait_associations"] + "/data.parquet",
        metadata_pq=config["trait_associations"] + "/metadata.parquet",
        metadata_tsv=config["trait_associations"] + "/metadata.tsv",
    input:
        phenotype_metadata_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
        genebass_pq=config["genebass_results_filtered_pq"],
        pred_score='/s/project/rep/processed/ukbb_wes_200k/pLoF_counts.parquet',
        plof_counts='/s/project/rep/processed/ukbb_wes_200k/pLoF_counts.parquet',
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