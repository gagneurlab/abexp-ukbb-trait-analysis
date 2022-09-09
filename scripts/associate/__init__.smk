SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

include: f"{SNAKEFILE_DIR}/associate_loftee_plof.py.smk"
include: f"{SNAKEFILE_DIR}/associate_per_tissue.py.smk"
include: f"{SNAKEFILE_DIR}/regression.py.smk"

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p


