SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''
ASSOCIATION_BASEPATH=f'''{COVARIATES_BASEPATH}/fset={{feature_set}}'''

OUTPUT_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}/fset={{feature_set}}/type={{model_type}}/polygenic_risk_score'''


rule associate__polygenic_risk_score:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        #prs_features_pq=directory(f"{OUTPUT_BASEPATH}/prs_features.parquet"),
        predictions_pq=f"{OUTPUT_BASEPATH}/predictions.parquet",
        r2_scores_pq=f"{OUTPUT_BASEPATH}/r2_scores.parquet",
        r2_scores_tsv=f"{OUTPUT_BASEPATH}/r2_scores.tsv",
        precision_recall_baseline_pq=f"{OUTPUT_BASEPATH}/precision_recall.baseline.parquet",
        precision_recall_full_pq=f"{OUTPUT_BASEPATH}/precision_recall.full.parquet",
        #prc_plot_png = f"{OUTPUT_BASEPATH}/precision_recall_plot.png",
        #predictions_plot_png = f"{OUTPUT_BASEPATH}/predictions_plot.png",
        # stats_pq=f"{OUTPUT_BASEPATH}/stats.parquet",
        # restricted_summary_txt=f"{OUTPUT_BASEPATH}/restricted_model.summary.txt",
        # full_summary_txt=f"{OUTPUT_BASEPATH}/full_model.summary.txt",
        # restricted_params_pq=f"{OUTPUT_BASEPATH}/restricted_model.params.parquet",
        # full_params_pq=f"{OUTPUT_BASEPATH}/full_model.params.parquet",
        # pickling not working bc of patsy formula
        # restricted_model=f"{OUTPUT_BASEPATH}/restricted_model.pickle",
        # full_model=f"{OUTPUT_BASEPATH}/full_model.pickle",
        # touch_file=touch(f"{OUTPUT_BASEPATH}/done"),
    input:
        # protein_coding_genes_pq=config["protein_coding_genes_pq"],
        associations_pq=f"{ASSOCIATION_BASEPATH}/associations.parquet",
        featureset_config=f"{ASSOCIATION_BASEPATH}/config.yaml",
        regression_done=f"{ASSOCIATION_BASEPATH}/done",
        # covariates
        covariates_pq=f'''{COVARIATES_BASEPATH}/covariates.parquet''',
        # clumping
        clumping_variants_pq=f'''{COVARIATES_BASEPATH}/clumping_variants.parquet''',
        # sample splits
        samples_pq=f"{COVARIATES_BASEPATH}/samples.parquet",
        samples_pq_done=f"{COVARIATES_BASEPATH}/samples.parquet.done",
        train_test_split_pq=f"{COVARIATES_BASEPATH}/train_test_split.parquet",
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        output_basedir=OUTPUT_BASEPATH,
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


del COVARIATES_BASEPATH
del ASSOCIATION_BASEPATH
del OUTPUT_BASEPATH

