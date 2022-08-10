SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


rule genes:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        protein_coding_genes_csv=config["protein_coding_genes_csv"],
        protein_coding_genes_pq=config["protein_coding_genes_pq"],
    input:
        gtf_file=config["gtf_file"],
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