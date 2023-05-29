#!/bin/bash

# install packages and check if other requirements are fulfilled
cat <<EOF | R --no-save --no-restore

if(!require(ggpubr)){
    install.packages("ggpubr", rev(.libPaths()))
    library(ggpubr)
}
if(!require(ukbtools)){
    install.packages("ukbtools", rev(.libPaths()))
    library(ukbtools)
}
library(arrow)
library(data.table)
library(stringr)

EOF
PIP_NO_DEPS=1 pip install "interpret==0.2.7" "interpret-core==0.2.7"