#!/bin/bash
#SBATCH -J 50Mpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=10:30:00
#SBATCH --export=ALL
#SBATCH -o ./log/50Mpc.o
#SBATCH --mail-type=end
#SBATCH --mail-user=jiaxuanl@princeton.edu
#SBATCH --account=jennyg

now=$(date +"%T")
echo "start time ... $now"
echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module load lsst

python simgal.py --distance=50 --age=1.0 --log_m_stars='[5.5, 6, 6.5, 7]'
python simgal.py --distance=50 --age=0.4 --log_m_stars='[5.5, 6, 6.5, 7]'

# python simgal.py --distance=20 --age=1.0 --log_m_stars='[5.0]'

# python simgal.py --distance=5 --age=1.0 --log_m_stars='[4.0]'

now=$(date +"%T")
echo "end time ... $now"
