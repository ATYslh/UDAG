#!/bin/bash
#SBATCH --job-name=sfcWind
#SBATCH --output=sfcWind.o%J
#SBATCH --error=sfcWind.o%J
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=philipp.heinrich@hereon.de
#SBATCH --mail-type=END,FAIL
#SBATCH --account=bb1203

srun python3 gather_data_nocmor_MEU3_EVAL.py