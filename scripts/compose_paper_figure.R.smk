SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


OUTPUT_BASEPATH=f'''{config["paper_figure"]}'''


rule compose_paper_figure:
    threads: 16
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        touch_file=touch(f"{OUTPUT_BASEPATH}/done"),
    input:
        **{
            **{
                f"compare_associations__{k}": expand(
                    v,
                    comparison="paper_figure",
                )
                for k, v in rules.compare_associations.output.items()
            },
            **{
                f"compare_risk_scores__{k}": expand(
                    v,
                    comparison="paper_figure",
                    model_type="lightgbm",
                )
                for k, v in rules.compare_risk_scores.output.items()
            },
            **{
                f"compare_risk_scores__{k}__linear": expand(
                    v,
                    comparison="paper_figure",
                    model_type="normalized_linear",
                )
                for k, v in rules.compare_risk_scores.output.items()
            },
            "Alanine_aminotransferase": rules.associate__plot_risk_score.output["plotting_done"].format(**{
                "phenotype_col": "Alanine_aminotransferase",
                "feature_set": "AbExp_all_tissues",
                "covariates": "sex_age_genPC_CLMP_PRS",
                "model_type": "lightgbm",
            }),
            "combined_qqplot_pq": rules.compare_associations.output["qq_plot_pq"].format(comparison="paper_figure_randomized"),
        }
    params:
        compare_associations_dir=rules.compare_associations.params.output_basedir.format(comparison="paper_figure"),
        compare_risk_scores_dir=rules.compare_risk_scores.params.output_basedir.format(
            comparison="paper_figure", 
            model_type="lightgbm",
        ),
        compare_risk_scores_linear_dir=rules.compare_risk_scores.params.output_basedir.format(
            comparison="paper_figure", 
            model_type="normalized_linear",
        ),
        alanine_aminotransferase_dir=rules.associate__plot_risk_score.params["output_basedir"].format(**{
            "phenotype_col": "Alanine_aminotransferase",
            "feature_set": "AbExp_all_tissues",
            "covariates": "sex_age_genPC_CLMP_PRS",
            "model_type": "lightgbm",
        }),
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        output_basedir=OUTPUT_BASEPATH,
    wildcard_constraints:
        covariates="[^/]+",
    conda: f'{CONDA_ENV_YAML_DIR}/ukbb-trait-analysis-R.yaml'
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.R"


del OUTPUT_BASEPATH

