#!/bin/bash
#SBATCH --account populationgenomics
#SBATCH --job-name=archaic_admixture
#SBATCH --output=archaic_admixture.out
#SBATCH --error=archaic_admixture.err
#SBATCH --time=00:30:00
#SBATCH -c 50
#SBATCH --mem-per-cpu=1GB

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hpc_course_r
Rscript archaic_admixture_EPAS1.R
