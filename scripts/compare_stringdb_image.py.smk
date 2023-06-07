SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


OUTPUT_BASEPATH=f'''{config["stringdb_comparison"]}/comp={{comparison}}/pheno={{phenotype}}:cov={{covariates}}'''
YAML_PATH=f"""{SNAKEFILE_DIR}/compare_stringdb_image@{{comparison}}.yaml"""


def _compare_stringdb_image_input_fn(wildcards, yaml_path=YAML_PATH):
    import yaml
    
    yaml_path = yaml_path.format(comparison=wildcards.comparison)
    
    with open(yaml_path, "r") as fd:
        config=yaml.safe_load(fd)
    
    input_files={
        **{
            k: expand(
                v,
                # feature_set_basepath + '/config.yaml',
                phenotype_col=wildcards["phenotype"], 
                feature_set=config["features_sets"],
                covariates=wildcards["covariates"],
            )
            for k, v in rules.associate__regression.output.items()
        },
        "most_associating_terms_pq": expand(
            rules.associate__compare_genebass.output["most_associating_terms_pq"],
            phenotype_col=wildcards["phenotype"], 
            feature_set=config["features_sets"],
            covariates=wildcards["covariates"],
        ),
        "significant_genes_pq": expand(
            rules.associate__compare_genebass.output["significant_genes_pq"],
            phenotype_col=wildcards["phenotype"], 
            feature_set=config["features_sets"],
            covariates=wildcards["covariates"],
        ),
    }

    return input_files



rule compare_stringdb_image:
    threads: 64
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        string_svg=f"{OUTPUT_BASEPATH}/stringdb_image.svg",
        string_png=f"{OUTPUT_BASEPATH}/stringdb_image.png",
        stringdb_ids_pq=f"{OUTPUT_BASEPATH}/stringdb_ids.parquet",
        stringdb_network_pq=f"{OUTPUT_BASEPATH}/stringdb_network.parquet",
        touch_file=touch(f"{OUTPUT_BASEPATH}/stringdb_image.done"),
    input: unpack(_compare_stringdb_image_input_fn)
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        output_basedir=OUTPUT_BASEPATH,
        config_yaml=YAML_PATH,
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
del YAML_PATH
