#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 10G
#SBATCH --time 1:00:00
#SBATCH --array=2-11
#SBATCH --job-name pre4Cin

gitHubDirectory=$1
FourCDirectory=$2
distanceToRemoveAroundVP=$3

# This script will smooth the segToFrag which are output of
# the 4Cseq analysis
path="$PWD/"
pathWithTableWithGenomes="${gitHubDirectory}/tables/table.txt"
pathWithBRFiles="${gitHubDirectory}/tables/"
pathForScripts="${gitHubDirectory}/scripts/"
pathForTable="${gitHubDirectory}/tables/4CinTable.txt"
pathWithResultsOf4C="${FourCDirectory}/analysisRodriguezCarballo2020/"

myChr="chr2"

module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0

# First get the info for the current experiment
if [ ! -e ${pathForTable} ]; then
  echo "${pathForTable} does not exists"
  exit 1
fi

genome=`cat ${pathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}'`
experiment=`cat ${pathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}'`
samples4CWithComma=`cat ${pathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $4}'`

if [ ! -e ${pathWithTableWithGenomes} ]; then
  echo "${pathWithTableWithGenomes} does not exists"
  exit 1
fi
brFiles=`cat ${pathWithTableWithGenomes} | awk -v g=${genome} '$1==g{print $2}'`
originalGenome=${genome}
if [ -z ${brFiles} ]; then
  assembly=${genome}
  if [ ${genome} = "mm10" ]; then
    genome="Wt"
  fi
else
  assembly=mm10_${genome}
fi

mkdir -p ${experiment}/raw_inputs
cd ${experiment}/raw_inputs

# Make a symlink of the segToFrag from the 4Cseq step in the raw_inputs folder
samples4C=`echo ${samples4CWithComma} | tr "," " "`
for sample in ${samples4C}; do
  if [ `ls ${pathWithResultsOf4C}/res_files_4Cseq_${genome}*/segToFrag_${sample}.bw | wc -l` -ne 1 ]; then
    echo "There is not exactly one segToFrag for ${sample} on ${genome}"
    exit 1
  fi
  ln -s ${pathWithResultsOf4C}/res_files_4Cseq_${genome}*/segToFrag_${sample}.bw .
done

# Prepare a config file for smoothing
template=${path}/${assembly}_NlaIIIDpnII.bedGraph.gz
pathForVPNameAndPos=${path}/${assembly}_viewpointPos_2col.txt
echo "### Required parameters
pathForFunctions4C <- \"${pathForScripts}/4C_functions.R\"
pathForVPNameAndPos <- \"${pathForVPNameAndPos}\"
myChr <- \"${myChr}\"
folderWithInputs <- \"$PWD\"

### Optional parameters
# Output location
outputFolder <- \"${path}${experiment}/smoothedBedGraphs_2kb_11frags\" #The output folder where will be the bedgraphs.

## Analysis parameters:
nbsOfFragmentPerWindow <- 11 #smoothing window

## Desired output:
outputFormat <- c(\"HTSStyle\")
# HTSstyle, it is asigned to the whole fragment: /!\ THIS MEANS IT IS NOT A \"NORMAL\" BEDGRAPH AND WE CANNOT USE THE CLASSICAL TOOLS TO EVALUATE THE SCORE IN A REGION.
smoothThrough <- T # Usually we do not want to spread signal from before the viewpoint after the viewpoint and vice-versa but it is necessary if you do not want to have a hole around the viewpoint.
distanceToRemoveAroundVP <- ${distanceToRemoveAroundVP} #Comment or put 0 to keep all information (for example to call peak)
smoothAfter <- T # By default, first you smooth the data, then you remove the signal around the viewpoint but for 4Cin you need to remove first and smooth after.
replaceNAby0<-T #IGV does not accept NA values. To make your bedGraphs compatible with IGV put T

### 4C PARAMETERS ###
template<-\"${template}\"
" > ${path}${experiment}/configFile_forSmooth_${experiment}.R

# If the template does not exists, create it
# if the index array is the first one on this genome:
# We need a template to know all the fragments possible.
# In the segToFrag, only the fragments with non 0 coverage are reported
# A template can be extracted from another output file of the 4Cseq step
# If you need such a template you can email me
firstIndexWithTheGenome=`cat ${pathForTable} | awk -v g=${originalGenome} '$1==g{print NR}' | head -n 1`
if [ "${SLURM_ARRAY_TASK_ID}" = "${firstIndexWithTheGenome}" ]; then
  if [ ! -e ${template} ]; then
    possibleT=`ls ${pathWithResultsOf4C}/res_files_4Cseq_${genome}*/normalised_scorePerFeature_*.bedGraph.gz | awk 'NR==1{print $1}'`
    if [ ! -z ${possibleT} ]; then
      zcat ${possibleT} | awk 'NF==4{print $1"\t"$2"\t"$3}' | gzip > ${template}
    else
      echo "There is no file ${genome}/normalised_scorePerFeature_*.bedGraph.gz"
      exit 1
    fi
  fi
  if [ ! -z ${brFiles} ]; then
    # This is not a WT genome so the file with coordinates may need to be created.
    if [ ! -e ${pathForVPNameAndPos} ]; then
      brFilesSpace=`echo ${brFiles} | tr "," " "`
      for f in ${brFilesSpace}; do
        if [ ! -e ${pathWithBRFiles}${f} ]; then
          echo "${pathWithBRFiles}${f} does not exists."
          exit 1
        fi
        ln -s ${pathWithBRFiles}${f} .
      done
      # We will convert the one of mm10 in bed and shift it
      cat ${path}/mm10_viewpointPos_2col.txt | awk -v mc=${myChr} '{print mc"\t"$2"\t"$2 + 1"\t"$1}' > temp.bed
      if [ ! -e ${pathForScripts}/shiftBedWithMultipleBR.sh ]; then
        echo "${pathForScripts}/shiftBedWithMultipleBR.sh does not exists"
        exit 1
      fi
      bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} temp.bed temp_${assembly}.bed 0
      cat temp_${assembly}.bed | awk '{print $4"\t"$2}' > ${pathForVPNameAndPos}
      rm temp.bed temp_${assembly}.bed
    fi
  fi
else
  # We wait that the one with the first index prepare the files.
  sleep 4m
fi

if [ ! -e ${template} ]; then
  echo "${template} does not exists but should have been created... this is strange"
  exit 1
fi

if [ ! -e ${pathForVPNameAndPos} ]; then
  echo "${pathForVPNameAndPos} does not exists but should have been created... this is strange"
  exit 1
fi

if [ ! -e ${pathForScripts}/smooth4CandUMI.R ]; then
  echo "${pathForScripts}/smooth4CandUMI.R does not exists"
  exit 1
fi

Rscript ${pathForScripts}/smooth4CandUMI.R ${path}${experiment}/configFile_forSmooth_${experiment}.R
