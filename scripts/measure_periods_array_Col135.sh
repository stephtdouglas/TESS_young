#!/bin/bash
#
#SBATCH --job-name=col135_prot
#SBATCH --output=/data/douglaslab/script_logs/slurm-%A_%a.out
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=2gb
#SBATCH --array=0-29
#SBATCH --requeue
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=douglste@lafayette.edu

srun hostname

source ~/.bashrc

cd ~/projects/TESS_young/
source activate lightkurve2
python measure_periods.py tables/Collinder_135_downloads_2021-06-17.csv Collinder_135

