library(ParallelStructure, lib.loc="/home/v1kdick3/R/x86_64-pc-linux-gnu-library/3.5/")
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run ParallelStructure                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

options(scipen=999)

# Genotypes

gen <- fread("/exports/eddie/scratch/v1kdick3/ADX_paper/4.analysis/Stru_1000/1.Input_files/ADX_geno95_TunNP.recode.strct_in_edit", header = F)

#~~ Specify in and out files for structure

infile <- "/exports/eddie/scratch/v1kdick3/ADX_paper/4.analysis/Stru_1000/1.Input_files/ADX_geno95_TunNP.recode.strct_in_edit"
outpath <- "/exports/eddie/scratch/v1kdick3/ADX_paper/4.analysis/Stru_1000/3.structure_out/"

#~~ construct job matrix and write to job file

nrep <- 10
up_to_k <- 10
niter <- 100000
burnin <- 50000

#~~ define variables for job matrix

k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

#~~ make the job matrix
pop <- "1,2,3,4,5,6,7"

jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
		                            rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(jobs), ncol = length(jobs[1,]), file = "/exports/eddie/scratch/v1kdick3/ADX_paper/4.analysis/Stru_1000/3.structure_out/jobs.txt")

#~~ file path to structure

STR_path='/exports/eddie/scratch/v1kdick3/tools/structure/'



#~~ Run structure

ParallelStructure::parallel_structure(structure_path=STR_path, 
				      joblist="/exports/eddie/scratch/v1kdick3/ADX_paper/4.analysis/Stru_1000/3.structure_out/jobs.txt", 
				      n_cpu=1,
				      infile=infile, 
				      outpath= outpath, 
#					  markernames = 0,
#					  mapdist = 0,
				      numinds = nrow(gen),
				      numloci= (ncol(gen)-2)/2,
				      locdata = 1,
				      popdata = 0, 
				      popflag = 0, 
				      usepopinfo = 0, 
				      locprior=0,
				      printqhat=1, 
				      plot_output=0, 
				      onerowperind=1)
