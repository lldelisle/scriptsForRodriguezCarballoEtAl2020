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

gitHubDirectory=$1

wd="$PWD/analysisRodriguezCarballo2020/"

module purge
module load sra-toolkit/2.9.6

mkdir -p ${wd}
cd ${wd}
# Download previously published demultiplexed fastq
sra=SRR5855209
fasterq-dump -o E12_PFL_wt_CS65_r1.fastq ${sra}

for ((i=1;i<=5;i++)); do
  sra="SRR$((5855184 + ${i}))"
  fasterq-dump -o E12_PFL_wt_CS38_r${i}.fastq ${sra}
done

for ((i=1;i<=3;i++)); do
  sra="SRR$((5855168 + ${i}))"
  fasterq-dump -o E12_PFL_wt_Hoxd4_r${i}.fastq ${sra}
done

# Download demultiplexed fastq of this analysis
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

# One is mapped on both
if [ -e invTDOM/E9_FLB_invTDOM_Hoxd11.fastq ]; then
  cp invTDOM/E9_FLB_invTDOM_Hoxd11.fastq Wt/
fi
