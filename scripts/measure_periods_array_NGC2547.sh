#!/bin/bash
#
#SBATCH --job-name=ngc2547_prot
#SBATCH --output=/data/douglaslab/script_logs/slurm-%A_%a.out
#SBATCH --qos=douglaslab
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=2gb
#SBATCH --array=0-17
#SBATCH --requeue
#
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=douglste@lafayette.edu

srun hostname

source ~/.bashrc

cd ~/projects/TESS_young/
source activate lightkurve2
python measure_periods.py tables/NGC_2547_downloads_2021-06-17.csv NGC_2547

