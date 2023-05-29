SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

include: f"{SNAKEFILE_DIR}/genes.py.smk"
include: f"{SNAKEFILE_DIR}/compare_associations.py.smk"
include: f"{SNAKEFILE_DIR}/compare_risk_scores.py.smk"
include: f"{SNAKEFILE_DIR}/compare_stringdb_image.py.smk"

# subdirectories
smkpaths = [
    # ukbb phenotypes:
    f"{SNAKEFILE_DIR}/ukbb_phenotypes/__init__.smk",
    # genebass:
    f"{SNAKEFILE_DIR}/genebass/__init__.smk",
    # feature sets for association tests
    f"{SNAKEFILE_DIR}/feature_sets/__init__.smk",
    # run association tests
    f"{SNAKEFILE_DIR}/associate/__init__.smk",
    # load the paper figure rule last
    f"{SNAKEFILE_DIR}/compose_paper_figure.R.smk"
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p


