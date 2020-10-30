#!/bin/bash -l

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 30G
#SBATCH --time 10:00:00
#SBATCH --job-name 4CinAve

gitHubDirectory=$1
plotDirectory=$2

# The current directory should be the parent folder of all 4Cin results
path="$PWD/"
pathForTable="${gitHubDirectory}/tables/4CinTable.txt"
pathForScripts="${gitHubDirectory}/scripts/"

module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0

# First make an archive with all vHiC:
tar -zcvmf allVhic.tar.gz 4Cin*/model*/wd/HoxD/*final*/vhic_HoxD.txt

# Then perform averages and clustering
Rscript ${pathForScripts}/correlations.R $plotDirectory

# We store the matrix of distance in a cool file
# We use a conda environment with cooler:
exists=`conda info --envs | awk '$1=="pgt_3.5"{print}' | wc -l`
if [ $exists -ne 1 ]; then
  conda create -y -n pgt_3.5 pygenometracks=3.5
fi
conda activate pgt_3.5

for f in average__*; do
  echo $f
  f_name=`basename $f .txt`
  experiment=`echo $f | awk '{split($1, a, "__"); print a[2]}'`
  dataFolder=`ls -d */model_${experiment}_onTDOM | head -n 1`
  python ${pathForScripts}/4CinToCool.py --data_dir ${dataFolder} --output_file ${f_name}.cool --pixel_file $f
done
