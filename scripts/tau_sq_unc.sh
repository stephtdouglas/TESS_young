#!/bin/bash
#
#SBATCH --job-name=tausq_syn
#SBATCH --output=/data/douglaslab/douglste/script_logs/slurm-%A-%a.out
#SBATCH --account=douglaslab
#SBATCH --partition=douglaslab,node
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:15:00
#SBATCH --mem-per-cpu=8gb
#SBATCH --array=20-59
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=douglste@lafayette.edu

srun hostname

source ~/.bashrc

cd /scratch/

mkdir -p douglste
cd douglste

mkdir -p plots

mkdir -p tables

cp ~/projects/TESS_young/tables/tausq_ZAMS_Compare*csv tables/
cp ~/projects/TESS_young/*py .
cp ~/projects/TESS_young/tab_*csv .

ulimit -n 8192

source activate lightkurve2

srun python tau_sq_synthetic_obs.py

mv tables/*SYN* ~/projects/TESS_young/tables/

mv plots/* ~/projects/TESS_young/plots/

cd ..

pwd

ls

ls ./*


