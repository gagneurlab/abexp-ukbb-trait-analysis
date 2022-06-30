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
if(!require(ukbtools)){
    install.packages("ukbtools", rev(.libPaths()))
    library(ukbtools)
}
library(arrow)
library(data.table)
library(stringr)
# library(ncdf4)

