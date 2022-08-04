#!/bin/bash
#
#SBATCH --job-name=ic2602_prot
#SBATCH --output=/data/douglaslab/script_logs/slurm-%A_%a.out
#SBATCH --account=douglaslab
#SBATCH --partition=douglaslab,node
#
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=2gb
#SBATCH --array=0-50%10
#SBATCH --requeue
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=douglste@lafayette.edu

srun hostname

source ~/.bashrc

cd ~/projects/TESS_young/
source activate lightkurve2
python measure_periods.py tables/IC_2602_downloads_2021-07-02.csv IC_2602

