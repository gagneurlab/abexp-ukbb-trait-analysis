SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

OUTPUT_BASEPATH=f'''{config["feature_sets_dir"]}/{SCRIPT}@agg={{agg}}'''

rule feature_sets__abexp:
    threads: 64
    resources:
        ntasks=1,
        mem_mb=250000
    output:
        data_pq=directory(f"{OUTPUT_BASEPATH}/data.parquet"),
        data_pq_done=touch(f"{OUTPUT_BASEPATH}/data.parquet.done"),
    input:
        abexp_predictions=config["abexp_predictions"],
        agg_config=ancient(f"{SNAKEFILE_DIR}/abexp@{{agg}}.yaml"),
    params:
        nb_script=f"{SCRIPT}",
#     log:
#         notebook=f"{DS_DIR}/feature_sets/{SCRIPT}@{{id}}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"

        
del OUTPUT_BASEPATH
