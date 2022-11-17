#!/bin/bash
#
#SBATCH --job-name=ic2391_prot
#SBATCH --output=/data/douglaslab/script_logs/slurm-%A_%a.out
#SBATCH --account=douglaslab
#SBATCH --partition=douglaslab,node
#
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=2gb
#SBATCH --array=0-15
#SBATCH --requeue
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=douglste@lafayette.edu

srun hostname

source ~/.bashrc

cd ~/projects/TESS_young/
source activate lightkurve2
python measure_periods.py tables/IC_2391_downloads_2021-06-21.csv IC_2391

