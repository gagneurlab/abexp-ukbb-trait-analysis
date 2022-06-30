# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: R [conda env:anaconda-drop2py]
#     language: R
#     name: conda-env-anaconda-drop2py-r
# ---

# %%
# if(!require(ukbtools)){
#     install.packages("ukbtools", .libPaths()[-1])
#     library(ukbtools)
# }
library(ukbtools)
library(arrow)
library(data.table)
library(stringr)
# library(ncdf4)


# %%
args = c(
    '/s/raw/ukbiobank/phenotypes/Phenotype_2018_ID22580',
    'ukb43546',
    '/s/project/uk_biobank/processed/phenotypes/Phenotype_2018_ID22580/ukb43546.data.parquet',
    '/s/project/uk_biobank/processed/phenotypes/Phenotype_2018_ID22580/ukb43546.meta.parquet'
)

# %%
args = commandArgs(trailingOnly = TRUE)

# %%
args

# %%
if (length(args) < 4) {
    stop("Usage: Rscript <this_script> <phenotype_dir> <ukbb_code> <output_data_pq> <output_metadata_pq>", call. = FALSE)
}

# %%
phenotype_dir = normalizePath(args[1])
ukbb_code = args[2]
output_data_pq = args[3]
output_metadata_pq = args[4]

# %%
# phenotype_dir = normalizePath("/s/raw/ukbiobank/phenotypes/Phenotype_2018_ID22580/")
# ukbb_code="ukb22580"
# output_data_pq = normalizePath("/s/project/uk_biobank/processed/phenotypes/ukb22580.data.parquet")
# output_metadata_pq = normalizePath("/s/project/uk_biobank/processed/phenotypes/ukb22580.meta.parquet")

# %%
# paths = Sys.glob('/s/raw/ukbiobank/phenotypes/*/*.enc_ukb')
# paths = str_match(paths, '(.*)/(.*)\\.enc_ukb')[,2:3]

# %%
ukbb_code

# %%
phenotype_dir

# %%
# To create a field code to name key
df_field <- ukb_df_field(ukbb_code, path = phenotype_dir)
df_field

# %%
write_parquet(df_field, file.path(output_metadata_pq))

# %%
# read table
df <- ukb_df(ukbb_code, path = phenotype_dir)
df

# %%
write_parquet(df, file.path(output_data_pq))


