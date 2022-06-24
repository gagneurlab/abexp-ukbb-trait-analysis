#'---
#' title: Significant genes per trait
#' author: Felix Brechtmann
#' wb:
#'  input:
#'  - genebass_associations: '/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/filtered_associations.tsv'
#'  output: '/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/signif_genes.tsv'
#' type: script
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---


library(data.table)

dt <- fread(snakemake@input[['genebass_associations']])
# dt <- fread('/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/filtered_associations.tsv')

bonferroni_cutoff <- 0.01 / nrow(dt[annotation == 'pLoF'])
bonferroni_cutoff


dt[annotation == 'pLoF' & Pvalue_Burden < bonferroni_cutoff, .N, by = description]

dt[annotation == 'pLoF' & description == 'HDL cholesterol' & Pvalue_Burden < bonferroni_cutoff]


dt[annotation == 'pLoF' & Pvalue_Burden < bonferroni_cutoff,
    .(gene_id, gene_symbol, description, phenocode)]

fwrite(dt[annotation == 'pLoF' & Pvalue_Burden < bonferroni_cutoff,
          .(gene_id, gene_symbol, description, phenocode)], 
       '/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/signif_genes.tsv', sep = '\t')