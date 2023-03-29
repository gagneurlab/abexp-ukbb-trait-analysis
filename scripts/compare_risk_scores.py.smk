SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


# OUTPUT_BASEPATH=f'''{config["trait_association_comparison"]}/cov={{covariates}}/fset={{feature_set}}'''
OUTPUT_BASEPATH=f'''{config["risk_score_comparison"]}/comp={{comparison}}'''
YAML_PATH=f"""{SNAKEFILE_DIR}/compare_risk_scores@{{comparison}}.yaml"""
COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''
FEATURE_SET_BASEPATH=f'''{COVARIATES_BASEPATH}/fset={{feature_set}}'''

def _compare_risk_scores_input_fn(wildcards, yaml_path=YAML_PATH, feature_set_basepath=FEATURE_SET_BASEPATH):
    import yaml
    
    yaml_path = yaml_path.format(comparison=wildcards.comparison)
    
    with open(yaml_path, "r") as fd:
        config=yaml.safe_load(fd)
    
    input_files={
        **{
            k: expand(
                v,
                # feature_set_basepath + '/config.yaml',
                phenotype_col=config["phenotypes"], 
                feature_set=config["features_sets"],
                covariates=config["covariates"],
            )
            for k, v in rules.associate__polygenic_risk_score.output.items()
        },
        "plot_risk_scores_done": expand(
            rules.associate__plot_risk_score.output["plotting_done"],
            phenotype_col=config["phenotypes"], 
            feature_set=config["features_sets"],
            covariates=config["covariates"],
        )
    }
    
    return input_files


rule compare_risk_scores:
    threads: 64
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        #significant_genes_pq=directory(f"{OUTPUT_BASEPATH}/significant_genes.parquet"),
        touch_file=touch(f"{OUTPUT_BASEPATH}/done"),
    input: unpack(_compare_risk_scores_input_fn)
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        output_basedir=OUTPUT_BASEPATH,
        config_yaml=YAML_PATH,
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


del OUTPUT_BASEPATH
del YAML_PATH
del FEATURE_SET_BASEPATH
del COVARIATES_BASEPATH
