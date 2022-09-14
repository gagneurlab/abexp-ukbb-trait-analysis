SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


COVARIATES_BASEPATH=f'''{config["trait_associations"]}/cov={{covariates}}'''
OUTPUT_BASEPATH=f'''{COVARIATES_BASEPATH}/fset={{feature_set}}'''
TEMPLATE_FILE=f"{SNAKEFILE_DIR}/regression@{{feature_set}}.yaml"


rule associate__regression:
    threads: 128
    resources:
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt
    output:
        associations_pq=directory(OUTPUT_BASEPATH + "/associations.parquet"),
        # config=f"{OUTPUT_BASEPATH}/config.yaml",
        touch_file=touch(f"{OUTPUT_BASEPATH}/done"),
    input:
        unpack(require_yaml_input(
            f"{OUTPUT_BASEPATH}/config.yaml",
        )), # additional input from featureset config yaml
        featureset_config=f"{OUTPUT_BASEPATH}/config.yaml",
        protein_coding_genes_pq=config["protein_coding_genes_pq"],
        # covariates
        covariates_pq=f'''{COVARIATES_BASEPATH}/covariates.parquet''',
        # clumping
        clumping_variants_pq=f'''{COVARIATES_BASEPATH}/clumping_variants.parquet''',
    params:
        nb_script=f"{SNAKEFILE_DIR}/{SCRIPT}",
        # age_col='age_when_attended_assessment_centre_f21003_0_0',
        # sex_col='sex_f31_0_0',
        # year_of_birth_col='year_of_birth_f34_0_0',
    wildcard_constraints:
        covariates="[^/]+",
#     log:
#         notebook=f"{DS_DIR}/{SCRIPT}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"


rule associate__regression_config:
    output:
        config=f"{OUTPUT_BASEPATH}/config.yaml"
    input:
        config_template=ancient(TEMPLATE_FILE),
        covariate_config=f'''{COVARIATES_BASEPATH}/config.yaml''',
    params:
        output_basedir=OUTPUT_BASEPATH,
#        age_col='age_when_attended_assessment_centre_f21003_0_0',
#        sex_col='sex_f31_0_0',
#        year_of_birth_col='year_of_birth_f34_0_0',
        nb_script=f"{SNAKEFILE_DIR}/regression_config.py",
    wildcard_constraints:
        ds_dir="[^/]+",
    script:
        "{params.nb_script}.py"
        

del COVARIATES_BASEPATH
del OUTPUT_BASEPATH
del TEMPLATE_FILE

localrules: associate__regression_config

