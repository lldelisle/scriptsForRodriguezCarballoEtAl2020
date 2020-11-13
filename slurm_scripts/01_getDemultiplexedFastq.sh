#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 10G
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --job-name getDemultiplexed

FourCAnalysis=$1
gitHubDirectory=$1

wd="$PWD/${FourCAnalysis}/"

module purge
module load sra-toolkit/2.9.6

mkdir -p ${wd}
cd ${wd}

# Download demultiplexed fastq
if [ ! -e ${gitHubDirectory}/tables/sraTable.txt ]; then
  echo "sraTable.txt does not exists"
  exit 1
fi
while read line; do
  sra=`echo ${line} | awk '{print $1}'`
  output=`echo ${line} | awk '{print $2".fastq"}'`
  fasterq-dump -o ${output} ${sra}
done < ${gitHubDirectory}/tables/sraTable.txt

# Put them in the folder accordingly to the genome on which to map
mkdir Wt
mkdir invTDOM

# All invTDOM and invTDOMdelBd are mapped on invTDOM
mv *invTDOM*.fastq invTDOM/
# All other are mapped on mm10
mv *.fastq Wt/

# All E9 invTDOM mapped on both
for f in invTDOM/E9*_invTDOM_*.fastq; do
  newName=`basename $f | awk '{gsub("invTDOM", "invTDOM_onwt", $1); print}'`
  cp $f Wt/${newName}
done
