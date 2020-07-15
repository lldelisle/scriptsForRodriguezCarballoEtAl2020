#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 24:00:00
#SBATCH --array=1-2
#SBATCH --job-name mutantGenome

gitHubDirectory=$1

path="$PWD/"
pathWithTableWithGenomes="${gitHubDirectory}/tables/table.txt"
pathForScripts="${gitHubDirectory}/scripts/"
pathWithBRFiles="${gitHubDirectory}/tables/"
firstEnzymeName="Nla"
secondEnzymeName="Dpn"
firstEnzymeRE="CATG"
secondEnzymeRE="GATC"
length=30

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load bowtie2/2.3.5
module load samtools/1.9
module load bedtools2/2.27.1
module load openblas/0.3.6-openmp 
module load r/3.6.0

if [ ! -e ${pathWithTableWithGenomes} ]; then
  echo "${pathWithTableWithGenomes} does not exists"
  exit 1
fi

genome=`cat ${pathWithTableWithGenomes} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}'`
brFiles=`cat ${pathWithTableWithGenomes} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}'`

if [ -z ${brFiles} ]; then
  assembly=${genome}
  if [ ${genome} = "mm10" ]; then
    genome="Wt"
  fi
else
  assembly=mm10_${genome}
fi

mkdir -p ${genome}

cd ${genome}

if [ -e ${assembly}.fa ]; then
  echo "genome already there"
else
  if [ ! -e ${path}/chr2_${genome}.fa ]; then
    if [ ${genome} = "Wt" ]; then
      echo "There is no chr2_Wt.fa in the path. You can download it from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/) make it all capital letters length 60 (for example with awk or seqtk seq -U for the upper case and fasta_formatter from fastx_toolkit or seqtk subseq for the 60bp length)."
    else
      echo "There is no chr2_${genome}.fa in the path. You can download it at https://doi.org/10.5281/zenodo.3826912"
    fi
    exit 1
  fi
  if [ ! -e ${path}/allChrsExceptchr2.fa.gz ]; then
    echo "There is no allChrsExceptchr2.fa.gz in the path. To make it, you need to download the sequence of chr1, chr3-chr19, chrX, chrY and chrM from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/), then concatenate them and make them all capital letters length 60 (for example with awk or seqtk seq -U for the upper case and fasta_formatter from fastx_toolkit or seqtk subseq for the 60bp length)."
    exit 1
  fi
  cp ${path}/allChrsExceptchr2.fa.gz ${assembly}.fa.gz
  gunzip ${assembly}.fa.gz
  cat ${path}/chr2_${genome}.fa >> ${assembly}.fa
fi

if [ ! -e ${assembly}.fa.fai ]; then
  samtools faidx ${assembly}.fa
fi

if [ ! -e ${assembly}.rev.1.bt2 ]; then
  bowtie2-build ${assembly}.fa ${assembly} &
fi

if [ ! -e ${assembly}.json ]; then
  if [ ! -e ${pathForScripts}/addNewGenomeToGenRep.R ]; then
    echo "${pathForScripts}/addNewGenomeToGenRep.R does not exists"
    exit 1
  fi
  Rscript ${pathForScripts}/addNewGenomeToGenRep.R ${assembly} ${assembly}.fa.fai ./
fi

if [ ! -e ${assembly}_rmsk.bed ]; then
  if [ ! -e ${path}mm10_rmsk.bed.gz ]; then
    echo "There is no ${path}mm10_rmsk.bed.gz. This file is obtained through UCSC website in tools table browser variations and repeat."
    exit 1
  fi
  if [ -z ${brFiles} ]; then
    gunzip -c ${path}mm10_rmsk.bed.gz > ${assembly}_rmsk.bed
  else
    brFilesSpace=`echo ${brFiles} | tr "," " "`
    for f in ${brFilesSpace}; do
      if [ ! -e ${pathWithBRFiles}${f} ]; then
        echo "${pathWithBRFiles}${f} does not exists but is required to convert the rmsk."
        exit 1
      fi
      ln -s ${pathWithBRFiles}${f} .
    done
    if [ ! -e ${pathForScripts}/shiftBedWithMultipleBR.sh ]; then
      echo "${pathForScripts}/shiftBedWithMultipleBR.sh does not exists"
      exit 1
    fi
    bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} ${path}mm10_rmsk.bed.gz ${assembly}_rmsk.bed 1
  fi
fi

libraryName=library_${genome}_${firstEnzymeName}${secondEnzymeName}_${length}bps

if [ ! -e ${libraryName}_segmentInfos.bed ]; then
  bedRptMaskFile=${assembly}_rmsk.bed
  if [ ! -e ${pathForScripts}/getRestEnzymeOccAndSeq.pl ]; then
    echo "${pathForScripts}/getRestEnzymeOccAndSeq.pl does not exists"
    exit 1
  fi
  if [ ! -e ${pathForScripts}/manual_createLibrary.sh ]; then
    echo "${pathForScripts}/manual_createLibrary.sh does not exists"
    exit 1
  fi
  ln -s ${pathForScripts}/getRestEnzymeOccAndSeq.pl .
  cat ${assembly}.fa | awk '{if($0~/^>/){print}else{print toupper($0)}}' | bash ${pathForScripts}/manual_createLibrary.sh -i - -m ${firstEnzymeRE} -s ${secondEnzymeRE} -l ${length} -r ${bedRptMaskFile} -n ${libraryName} > ${libraryName}.log 2>${libraryName}.err
fi
wait
