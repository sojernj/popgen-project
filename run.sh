#!/bin/bash
#SBATCH --account populationgenomics
#SBATCH --job-name=archaic_admixture
#SBATCH --output=archaic_admixture.out
#SBATCH --error=archaic_admixture.err
#SBATCH --time=02:00:00
#SBATCH -c 70
#SBATCH --mem-per-cpu=2GB

source $(conda info --base)/etc/profile.d/conda.sh
conda activate hpc_course_r
Rscript archaic_admixture.R