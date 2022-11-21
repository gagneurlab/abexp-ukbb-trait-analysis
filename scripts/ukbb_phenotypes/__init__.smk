SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]


# ukbb phenotypes
include: f"{SNAKEFILE_DIR}/read_phenotypes.R.smk"
include: f"{SNAKEFILE_DIR}/merge_phenotype_metadata.py.smk"
include: f"{SNAKEFILE_DIR}/merge_all_phenotypes.py.smk"
include: f"{SNAKEFILE_DIR}/decode_phenotype.py.smk"
