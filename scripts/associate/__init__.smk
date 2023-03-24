SNAKEFILE = workflow.included_stack[-1]
SNAKEFILE_DIR = os.path.dirname(SNAKEFILE)

SCRIPT=os.path.basename(SNAKEFILE)[:-4]

include: f"{SNAKEFILE_DIR}/covariates.py.smk"
# include: f"{SNAKEFILE_DIR}/associate_per_tissue.py.smk"
include: f"{SNAKEFILE_DIR}/regression.py.smk"
include: f"{SNAKEFILE_DIR}/compare_genebass.py.smk"
include: f"{SNAKEFILE_DIR}/compare_params.py.smk"
include: f"{SNAKEFILE_DIR}/qq_plot.py.smk"
include: f"{SNAKEFILE_DIR}/polygenic_risk_score.py.smk"
include: f"{SNAKEFILE_DIR}/plot_risk_score.py.smk"
# include: f"{SNAKEFILE_DIR}/lm_plot.py.smk"
include: f"{SNAKEFILE_DIR}/samples.py.smk"
include: f"{SNAKEFILE_DIR}/train_test_split.py.smk"

# subdirectories
smkpaths = [
]

for p in smkpaths:
    eprint("Including '%s'..." % p)
    include: p


