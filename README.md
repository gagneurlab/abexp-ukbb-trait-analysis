# UK Biobank trait analysis for AbExp

dev setup:
1) clone this repo and `cd` into it
2) run `find scripts/ -iname "*[.py.py|.R.R]" -exec conda run --name jupyterhub jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks
3) execute `run_slurm_jobs.sh --rerun-incomplete --rerun-triggers mtime -k --restart-times 1` to run the snakemake pipeline
