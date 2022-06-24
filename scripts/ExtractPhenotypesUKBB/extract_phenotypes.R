
# read traits
library(data.table)
traits <- fread('/s/project/bayesRare/UKBB_Splicing_Analysis/input/traits_monti.tsv')

traits[, 1:3]
phenocodes <- traits[, PHENOCODE]

# additionally retrive covariates age 21003, sex 31
phenocodes <- c(phenocodes, 21003, 31)


get_columns <- function(phenotypes){
    selection <- c()
    for(p in phenocodes){
      selection <- c(selection, grep(paste0("f.", p, '.0.0'), colnames(phenotypes)))
    }
    cols <- c("f.eid", colnames(phenotypes)[selection])
    
    phenotypes_trait_sub <- phenotypes[, ..cols]
    
    colnames(phenotypes_trait_sub) <- sapply(strsplit(colnames(phenotypes_trait_sub), split = '\\.'), function(x){ x[[2]]})
    phenotypes_trait_sub
}


# read latest phenotype file.
phenotypes1 <- fread('/s/raw/ukbiobank/phenotypes/Phenotype_20210327_ID46901/ukb46901.tab')
phenotypes1 <- get_columns(phenotypes1)

fwrite(phenotypes1, '/s/project/bayesRare/UKBB_Splicing_Analysis/ukbb_phenotypes.tsv', sep = '\t') 
