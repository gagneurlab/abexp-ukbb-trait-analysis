SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

OUTPUT_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}/fset={{feature_set}}'''


rule associate__get_stringdb_image:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        string_svg=f"{OUTPUT_BASEPATH}/stringdb_image.svg",
        string_png=f"{OUTPUT_BASEPATH}/stringdb_image.png",
        stringdb_ids_pq=f"{OUTPUT_BASEPATH}/stringdb_ids.parquet",
        stringdb_network_pq=f"{OUTPUT_BASEPATH}/stringdb_network.parquet",
        stringdb_enrichment_pq=f"{OUTPUT_BASEPATH}/stringdb_enrichment.parquet",
        touch_file=touch(f"{OUTPUT_BASEPATH}/stringdb_image.done"),
    input:
        significant_genes_pq=f"{OUTPUT_BASEPATH}/significant_genes.parquet",
        # significant_genes_tsv=f"{OUTPUT_BASEPATH}/significant_genes.tsv",
        # newly_found_genes_tsv=f"{OUTPUT_BASEPATH}/newly_found_genes.tsv",
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        output_basedir=OUTPUT_BASEPATH,
        pval_cutoff=0.05,
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


del OUTPUT_BASEPATH


