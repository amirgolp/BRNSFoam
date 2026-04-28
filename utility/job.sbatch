#!/bin/bash -l
#SBATCH --job-name=lowFlux
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --output=BRNS_%j.out
#SBATCH --error=BRNS_%j.err
#SBATCH --constraint=[CPU_ARCH:avx2|CPU_ARCH:avx512]

cd "$SLURM_SUBMIT_DIR"

module purge
module load       GCC/11.2.0  OpenMPI/4.1.1
module load OpenFOAM/v2212
source $FOAM_BASH
source /home/nhdfamgo/works/GeoChemFoam-5.1/etc/bashrc

# srun --cpu_bind=cores --distribution=block:cyclic interBRNSFoam
# mpiexec -np 6 interBRNSFoam -parallel
interBRNSFoam

reconstructPar
