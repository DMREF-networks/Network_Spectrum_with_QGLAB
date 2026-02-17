#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=64g
#SBATCH -t 0-12:59:59

module purge
module add python/3.12.4
python NetGen.py $1 $2 $3 $4 $5 $7 $8 

module purge
module add matlab/2025a

#matlab -singleCompThread -nosplash -nodesktop -r "open('Quantum-Graphs-master/QGObject.prj')"
cd Quantum-Graphs-master
matlab -singleCompThread -nosplash -nodesktop -r "SpectrumScript $1 $2 $3 $4 $5 $6 $7 $8"
