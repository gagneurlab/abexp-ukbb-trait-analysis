SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


OUTPUT_BASEPATH=f'''{config["trait_associations"]}/{SCRIPT}'''


rule associate_per_tissue:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        associations_pq=OUTPUT_BASEPATH + "/associations.parquet",
        data_pq=OUTPUT_BASEPATH + "/data.parquet",
        metadata_pq=OUTPUT_BASEPATH + "/metadata.parquet",
        metadata_tsv=OUTPUT_BASEPATH + "/metadata.tsv",
    input:
        phenotype_metadata_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
        genebass_pq=config["genebass_results_filtered_pq"],
        abexp_predictions=config["abexp_predictions"],
        plof_counts='/s/project/rep/processed/ukbb_wes_200k/pLoF_counts.parquet',
        protein_coding_genes_pq=config["protein_coding_genes_pq"],
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

del OUTPUT_BASEPATH