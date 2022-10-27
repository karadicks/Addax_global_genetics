#!/bin/sh
# Grid Engine options
#$ -N parallel_stru
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -R y
#$ -pe sharedmem 10

# Jobscript to run Parallel Structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load R

module load R

# Run

R CMD BATCH run_structure.R 


