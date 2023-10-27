#!/bin/bash

# GWAS imputation prerequisites 
# Emilio Ugalde
# Jan 27, 2022
# ssh emiuga@vector.meb.ki.se
# for guiding: https://packaging.python.org/en/latest/tutorials/installing-packages/

# PREREQUISITES
# https://github.com/hakyimlab/summary-gwas-imputation
# The basic requirements for running GWAS summary-imputation are python>=3.5 with the following packages:

# working directory
cd /nfs/GENETEC/TWAS/

# Clone repository
DIR=Programs/summary-gwas-imputation
mkdir -p $DIR
git clone https://github.com/hakyimlab/summary-gwas-imputation.git $DIR

# set up virtual environment 
python3 -m venv $DIR
source $DIR/bin/activate

# Intall libraries
pip3 install pandas==0.25.3
pip3 install scipy==1.4.1
pip3 install numpy==1.18.1
pip3 install bgen_reader==3.0.2
pip3 install cyvcf2==0.20.0
pip3 install pyliftover==0.4
pip3 install pyarrow==0.11.1



