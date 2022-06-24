#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import hail as hl

traits = pd.read_csv(snakemake.input['traits'], sep = "\t")
mt = hl.read_matrix_table(snakemake.input['genebass'])

mt.head(n=5)
traits = hl.literal(traits.PHENOCODE.astype('str').to_list())
traits

mt.filter_cols(traits.contains(mt.phenocode)
               ).entries().drop('interval',
                                'markerIDs',
                                'markerAFs',
                                "category").export(snakemake.output['filtered_genebass'])

