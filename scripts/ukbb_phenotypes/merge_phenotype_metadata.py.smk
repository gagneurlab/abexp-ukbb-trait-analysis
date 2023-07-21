SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

if not UKBB_SKIP_METADATA_UPDATE:
    rule merge_phenotype_metadata:
        threads: 16
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
        output:
            merged_meta_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/merged.meta.parquet",
            latest_meta_pq=f"{UKBB_PROCESSED_PHENOTYPES_DIR}/latest.meta.parquet",
        input:
            pheno_meta_pqs = expand(UKBB_PROCESSED_PHENOTYPES_DIR + "/{phenotype_dir}/{ukbb_code}.meta.parquet", zip, phenotype_dir=phenotype_dirs, ukbb_code=ukbb_codes),
        params:
            nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
            path_codes=list(zip(
                expand(UKBB_PROCESSED_PHENOTYPES_DIR + "/{phenotype_dir}/{ukbb_code}.meta.parquet", zip, phenotype_dir=phenotype_dirs, ukbb_code=ukbb_codes),
                phenotype_dirs,
                ukbb_codes,
            ))
        wildcard_constraints:
            ds_dir="[^/]+",
    #     log:
    #         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
    #     notebook:
    #         "{params.nb_script}.ipynb"
        script:
            "{params.nb_script}.py"
