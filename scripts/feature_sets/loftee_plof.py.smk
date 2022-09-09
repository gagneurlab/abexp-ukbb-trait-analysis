SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

import yaml

OUTPUT_BASEPATH=f'''{config["feature_sets_dir"]}/{SCRIPT}'''

rule feature_sets__loftee_plof:
    threads: 64
    resources:
        ntasks=1,
        mem_mb=250000
    output:
        data_pq=directory(f"{OUTPUT_BASEPATH}/data.parquet"),
        data_pq_done=touch(f"{OUTPUT_BASEPATH}/data.parquet.done"),
    input:
        vep_predictions=config["vep_predictions"],
    params:
        nb_script=f"{SCRIPT}",
#     log:
#         notebook=f"{DS_DIR}/feature_sets/{SCRIPT}@{{id}}.ipynb"
#     notebook:
#         "{params.nb_script}.ipynb"
    script:
        "{params.nb_script}.py"

        
del OUTPUT_BASEPATH
