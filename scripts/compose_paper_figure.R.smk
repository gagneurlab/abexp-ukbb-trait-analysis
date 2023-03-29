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
                )
                for k, v in rules.compare_risk_scores.output.items()
            },
            "Alanine_aminotransferase": rules.associate__plot_risk_score.output["plotting_done"].format(**{
                "phenotype_col": "Alanine_aminotransferase",
                "feature_set": "AbExp_all_tissues",
                "covariates": "sex_age_genPC_CLMP_PRS"
            }),
        }
    params:
        compare_associations_dir=rules.compare_associations.params.output_basedir.format(comparison="paper_figure"),
        compare_risk_scores_dir=rules.compare_risk_scores.params.output_basedir.format(comparison="paper_figure"),
        alanine_aminotransferase_dir=rules.associate__plot_risk_score.params["output_basedir"].format(**{
            "phenotype_col": "Alanine_aminotransferase",
            "feature_set": "AbExp_all_tissues",
            "covariates": "sex_age_genPC_CLMP_PRS"
        }),
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


del OUTPUT_BASEPATH

