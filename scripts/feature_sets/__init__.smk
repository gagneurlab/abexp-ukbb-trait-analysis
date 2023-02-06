SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

include: f"{SNAKEFILE_DIR}/loftee_plof.py.smk"
include: f"{SNAKEFILE_DIR}/abexp.py.smk"
include: f"{SNAKEFILE_DIR}/abexp_pivot.py.smk"
include: f"{SNAKEFILE_DIR}/merged_abexp_pivot_loftee.py.smk"

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p


