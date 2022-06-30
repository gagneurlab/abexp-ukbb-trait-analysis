SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

include: f"{SNAKEFILE_DIR}/associate.py.smk"

# subdirectories
smkpaths = [
    # ukbb phenotypes:
    f"{SNAKEFILE_DIR}/ukbb_phenotypes/__init__.smk",
    # genebass:
    f"{SNAKEFILE_DIR}/genebass/__init__.smk",
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p


