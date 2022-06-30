SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


# genebass
include: f"{SNAKEFILE_DIR}/extract_genebass.smk"
include: f"{SNAKEFILE_DIR}/filter_genebass.py.smk"
