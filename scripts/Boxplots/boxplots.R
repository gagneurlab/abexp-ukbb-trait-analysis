#'---
#' title: Boxplot per trait
#' author: Felix Brechtmann
#' wb:
#'  input:
#'  - genebass_associations: '/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/filtered_associations.tsv'
#'  output: 
#'  - plot: '/s/project/bayesRare/UKBB_Splicing_Analysis/results/{phenocode}/boxplot.png'
#'  - wBhtml: '`sm config["htmlOutputPath"] + "/results/{phenocode}/boxplot.html"`'
#'  type: noindex
#' type: script
#' output:
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---




library(ggplot2)
library(data.table)
library(ggthemes)
library(cowplot)

phenoc <- snakemake@wildcards[['phenocode']]
#phenoc <- 30760
phenocColumn <- paste0('f', phenoc)

# read table of associated gene
associated_genes <- fread('/s/project/bayesRare/UKBB_Splicing_Analysis/genebass/signif_genes.tsv')
# subset to phenocode
associated_genes <- associated_genes[phenocode == phenoc]

# read phenotype info
pheno <- fread('/s/project/bayesRare/UKBB_Splicing_Analysis/ukbb_phenotypes.tsv')
colnames(pheno)[2:length(colnames(pheno))] <- paste0('f',colnames(pheno)[2:length(colnames(pheno))])
colnames(pheno)

#keep sex, age, sample id and phenotype
col_of_interest <- c('eid', 'f21003', 'f31',  paste0('f', phenoc))
pheno <- pheno[, ..col_of_interest]

# remove NA values from the phenotype column
pheno <- pheno[!is.na(get(phenocColumn))]

# regress out standard covariates age (f21003), age ^ 2, sex (f31) , age × sex
m <- lm(paste0(phenocColumn, '~ f31 + f21003 + I(f21003 ** 2) + f31 * f21003'), data = pheno)
pheno[, residuals := residuals(m)]

ggplot(pheno, aes(residuals)) + geom_histogram()

# read number of pLoF per gene
pLoF <- fread('/s/project/bayesRare/rarephenopred/misc/pLoF_counts_multitrait.tsv')

meltpLoF <- melt(pLoF, id.vars = 'sample_id', variable.name = 'gene_id', value.name = 'pLoF')
meltpLoF <- meltpLoF[gene_id %in% associated_genes[, gene_id]]

# read in splicing variants
splice_vars <- fread('/s/project/bayesRare/rarephenopred/misc/splicing_variants_multitrait.tsv')
splice_vars <- splice_vars[gene_id %in% associated_genes[, gene_id]]

splice_vars_gene_Splice_AI <- splice_vars[delta_score >= 0.8, .(SpliceAI = sum(GT)), by=c('gene_id', 'sample_id') ]
splice_vars_gene_Ab_Splice <- splice_vars[AbSplice_DNA__max >= 0.075, .(AbSplice = sum(GT)), by=c('gene_id', 'sample_id') ]
#splice_vars_gene_Ab_Splice <- splice_vars[gtex-grch37-adipose_subcutaneous.AbSplice_DNA > 0.075, .(AbSplice = sum(GT)), by=c('gene_id', 'sample_id') ]

splice_vars_gene <- merge(splice_vars_gene_Splice_AI, splice_vars_gene_Ab_Splice, all = TRUE)
splice_vars_gene[, table(SpliceAI, AbSplice)]

all_vars <- merge(meltpLoF, splice_vars_gene, by = c('sample_id', 'gene_id'), all = TRUE)

all_var_counts <- melt(all_vars, measure.vars = c('pLoF','SpliceAI', 'AbSplice'), 
                       variable.name = 'VariantType', 
                       value.name = 'Num_var')

all_var_counts[is.na(Num_var), Num_var := 0]


all_var_counts[, Num_vars_in_sample := sum(Num_var), by = sample_id]
all_var_counts[ Num_vars_in_sample == 0, uniqueN(sample_id)]
all_var_counts[, uniqueN(sample_id)]


all_var_counts[, num_pLoF := sum(VariantType == 'pLoF' & Num_var > 0), by = gene_id]
all_var_counts[, num_AbSplice := sum(VariantType == 'AbSplice' & Num_var > 0), by = gene_id]
all_var_counts[, num_spliceAI := sum(VariantType == 'SpliceAI' & Num_var > 0), by = gene_id]
all_var_counts[, Gene := paste0(gene_id, "\n N pLoF variants = ", num_pLoF,
                          "\n N SpliceAI variants = ",num_spliceAI,
                          "\n N AbSplice variants = ",num_AbSplice
)]
all_var_counts[, Rare_Variant := Num_var > 0]
all_var_counts

# merge with phenotype 
plot_table <- merge(all_var_counts, pheno[!is.na(residuals)], by.x = 'sample_id', by.y = 'eid')



boxplot <- ggplot(plot_table, aes(Gene, residuals, color = VariantType, linetype = Rare_Variant)) + geom_boxplot() +
  theme_cowplot() + background_grid() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
                                              legend.position = 'bottom') 

print(associated_genes[, unique(description)])

#+ fig.height=12, fig.width=8

boxplot


ggsave(snakemake@output[['plot']], boxplot)