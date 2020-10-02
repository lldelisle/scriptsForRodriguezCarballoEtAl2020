#!/bin/bash -l

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 30G
#SBATCH --time 48:00:00
#SBATCH --job-name 4Cin
#SBATCH --array=2-11

gitHubDirectory=$1

path="$PWD/"
pathForTable="${gitHubDirectory}/tables/4CinTable.txt"
pathWithTableWithGenomes="${gitHubDirectory}/tables/table.txt"
pathWithBRFiles="${gitHubDirectory}/tables/"
pathForScripts="${gitHubDirectory}/scripts/"

myChr="chr2"
start=74400000
end=75800000

module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0

# The current directory should ends by a letter but the inputs are in the folder without the letter
pathWithInputs=`echo $PWD | awk '{print substr($1, 1, length($1) - 1)}'`
cp -r ${pathWithInputs}/* .

if [ ! -e ${pathForTable} ]; then
  echo "${pathForTable} does not exists"
  exit 1
fi

genome=`cat ${pathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $1}'`
experiment=`cat ${pathForTable} | awk -v i=${SLURM_ARRAY_TASK_ID} 'NR==i{print $2}'`

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

pathForVPandPos="${path}/${assembly}_viewpointPos_2col.txt"

smallName=${experiment}_onTDOM

inputFolder="${path}${experiment}/smoothedBedGraphs_2kb_11frags/"
dataFolder="${path}/model_${smallName}/"
workingDir="${dataFolder}/wd/"

if [ -e ${dataFolder} ]; then
  rm -r ${dataFolder}
fi


if [ ! -z ${brFiles} ]; then
  # This is not a WT genome the region need to be adjusted
  # To avoid similar jobs to collapse
  # We will use a temporary folder
  temp_folder=`mktemp -d`
  cd ${temp_folder}
  brFilesSpace=`echo ${brFiles} | tr "," " "`
  for f in ${brFilesSpace}; do
    if [ ! -e ${pathWithBRFiles}${f} ]; then
      echo "${pathWithBRFiles}${f} does not exists."
      exit 1
    fi
    ln -s ${pathWithBRFiles}${f} .
  done
  # The region is converted to bed
  echo -e "${myChr}\t${start}\t${end}" > temp.bed
  # And is shifted
  if [ ! -e ${pathForScripts}/shiftBedWithMultipleBR.sh ]; then
    echo "${pathForScripts}/shiftBedWithMultipleBR.sh does not exists"
    exit 1
  fi
  bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} temp.bed temp_${genome}.bed 0
  # In case of multiple intervals the start is the minimum of all intervals
  start=`cat temp_${genome}.bed | awk 'BEGIN{min=1000000000}{if($2<min){min=$2}}END{print min}'`
  # end is maximum
  end=`cat temp_${genome}.bed | awk '{if($3>max){max=$3}}END{print max}'`
  # Go back to path
  cd $path
fi


mkdir -p ${workingDir}

inputFolderSmall=${inputFolder}/${start}_${end}
if [ -e ${inputFolderSmall} ]; then
  rm -r ${inputFolderSmall}
fi
mkdir -p ${inputFolderSmall}

# The bedgraph are restricted to the start end
if [ `ls ${inputFolder}/*.bedGraph | wc -l` -eq 0 ]; then
  echo "There is no bedGraph in ${inputFolder}"
  exit 1
fi
for f in ${inputFolder}/*.bedGraph; do
  name=`basename ${f} | awk '{split($0,a,"_smoothed");print a[1]}'`
  cat ${f} | awk -v s=${start} -v e=${end} '$2<e&&$3>s{print}' > ${inputFolderSmall}/${name}
done

# We run the script prepare_data.py from 4Cin
if [ ! -e ${pathForScripts}/prepare_data.py ]; then
  echo "${pathForScripts}/prepare_data.py does not exists"
  exit 1
fi
python ${pathForScripts}/prepare_data.py ${inputFolderSmall}/*

mkdir -p ${dataFolder}

for f in ${inputFolderSmall}/*_modified; do
  name=`basename ${f} | awk '{split($0,a,"_modified");print a[1]}'`
  # We assume that the name contains the name of the viewpoint separated by _
  # We will get the line
  line=`cat ${pathForVPandPos} | awk -v n=${name} '
    BEGIN{
      split(n, ns, "_")
      for (i in ns) keywords[ns[i]] = ""
    }
    $1 in keywords{
      print $0
    }'`
  if [ -z "${line}" ]; then
    echo "No corresponding line was found for ${name} in ${pathForVPandPos}"
    exit 1
  fi
  vp=`echo ${line} | awk '{print $1}'`
  pos=`echo ${line} | awk '{print $2}'`
  mv ${f} ${dataFolder}/${vp}
  echo "${vp} ${myChr}:${pos}" >> ${dataFolder}/primers.txt
done

cd ${workingDir}

# We use the docker image to run 4Cin
shifter --image batxes_fromdockerhub/4cin_ubuntu 4Cin.py --working_dir ${workingDir} --cpu 28 --verbose ${dataFolder} HoxD 

# Sometimes it stops at the preprocessing step...
if [ `ls HoxD/ | grep HoxD_output | wc -l` -eq 0 ]; then
  # You need to clean
  rm -r pre_modeling
  # You need to relaunch it
  shifter --image batxes_fromdockerhub/4cin_ubuntu 4Cin.py --working_dir ${workingDir} --cpu 28 --verbose ${dataFolder} HoxD 
fi
# Sometimes it cannot go to the end because of way to print float...
# For example 4Cin will raise:
# OSError: [Errno 2] No such file or directory: '/scratch/ldelisle/all4Cin/model_E9_FLB_wt_like_invCS3840_onTDOM/wd/HoxD/HoxD_output_0.3_-0.1_12000/'
# Whereas you have HoxD_output_0.30000000000000004_-0.1_12000
uZ=`ls HoxD/ | grep HoxD_output | awk -F "_" '{printf("%.1f",$(NF-2))}'`
lZ=`ls HoxD/ | grep HoxD_output | awk -F "_" '{printf("%.1f",$(NF-1))}'`
maxD=`ls HoxD/ | grep HoxD_output | awk -F "_" '{print $NF}'`

if [ ! -e HoxD/HoxD_final_output* ]; then
  # There was an error of float print
  mv HoxD/HoxD_output* HoxD/HoxD_output_${uZ}_${lZ}_$maxD
  shifter --image batxes_fromdockerhub/4cin_ubuntu 4Cin.py --working_dir $workingDir --cpu 28 --verbose $dataFolder HoxD --jump_step modeling --max_distance $maxD --uZ $uZ --lZ $lZ
fi

if [ ! -e HoD/log.txt ]; then
  shifter --image batxes_fromdockerhub/4cin_ubuntu 4Cin.py --working_dir $workingDir --verbose $dataFolder HoxD --jump_step analysis --max_distance $maxD --uZ $uZ --lZ $lZ
fi

cd ${path}
# We store the matrix of distance in a cool file
# We use a conda environment with cooler:
exists=`conda info --envs | awk '$1=="pgt_3.5"{print}' | wc -l`
if [ $exists -ne 1 ]; then
  conda env create -n pgt_3.5 pygenometracks=3.5  
fi
conda activate pgt_3.5

if [ ! -e ${pathForScripts}/4CinToCool.py ]; then
  echo "${pathForScripts}/4CinToCool.py does not exists"
  exit 1
fi
python ${pathForScripts}/4CinToCool.py --data_dir ${dataFolder} --output_file $(basename ${dataFolder}).cool
