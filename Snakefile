import pathlib
import wbuild
config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

configfile: "wbuild.yaml"
include: config['wBuildPath'] + "/wBuild.snakefile"
include: 'additional_rules/rules/data_preprocessing.smk'

import pandas as pd
associations = pd.read_csv('/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/signif_genes.tsv', sep = '\t')

rule all:
	input: rules.Index.output, 
	       expand('/s/project/bayesRare/UKBB_Splicing_Analysis/results/{phenocode}/boxplot.png', phenocode = set(associations.phenocode))
