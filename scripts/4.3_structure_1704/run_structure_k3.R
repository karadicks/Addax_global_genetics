#! /usr/bin/Rscript --vanilla --default-packages=utils

library(ParallelStructure)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run ParallelStructure                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



options(scipen=999)

setwd("/mnt/profiles/kdicks/3.WG_Projects/WG1705_ADX/2021_ADX_paper_analysis/4.analysis/3.structure_1704/")

#~~ Specify in and out files for structure

job="3"



infile <- "1.input_files/ADX_geno95_TunNP.recode.strct_in"
outpath <- "2.output_files/"

# Genotypes

gen <- fread(infile, header = F)

#~~ file path to structure

STR_path='/mnt/profiles/kdicks/tools/structure/'


#~~ Run Parallel Structure

ParallelStructure::parallel_structure(structure_path=STR_path, 
                                      joblist=paste0("2.output_files/stru_jobs_",job,".txt"), 
                                      n_cpu=15,
                                      infile=infile, 
                                      outpath= outpath, 
                                      numinds = nrow(gen),
                                      numloci=(ncol(gen)-2)/2, 
                                      popdata = 1, 
                                      #popflag = 1, 
                                      #usepopinfo = 1, 
                                      printqhat=1, 
                                      plot_output=0, 
                                      onerowperind=1)



