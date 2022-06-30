#!/usr/bin/env python
# coding: utf-8

import os

import pandas as pd
import hail as hl

### --- setup spark/hail --- ###
MEM = os.popen("ulimit -m").read()
if MEM.startswith("unlimited"):
    print("Memory not constrained, using all available memory...")
    import psutil
    MEM = psutil.virtual_memory().available / 1024

# reduce driver memory to 80% of available memory and convert to integer
MEM = int(MEM * 0.8)

print("memory: %dk" % MEM)

os.environ['PYSPARK_SUBMIT_ARGS'] = " ".join([
    '--driver-memory %dk' % MEM,
    'pyspark-shell'
])
os.environ['PYSPARK_SUBMIT_ARGS']

# init hail
hl.init(tmp_dir=os.environ.get("TMP"))

### --- setup finished --- ###

RESULTS_MT = snakemake.input["results_mt"]
RESULTS_PQ = snakemake.output["results_pq"]

mt = hl.read_matrix_table(
    RESULTS_MT
#     snakemake.input['genebass']
)


# fetch gene entries
mt_entries = (
    mt
    .entries()
    .drop(
        'interval',
        'markerIDs',
        'markerAFs',
        "category",
    )
)

# export as parquet
(
    mt_entries
    .to_spark()
    .write
    .parquet(RESULTS_PQ, mode="overwrite")
)


# read traits of interest
# traits = pd.read_csv(snakemake.input['traits'], sep = "\t")
# traits = hl.literal(traits.PHENOCODE.astype('str').to_list())
# traits