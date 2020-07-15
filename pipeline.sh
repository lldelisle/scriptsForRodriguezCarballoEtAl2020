# Define all paths
gitHubDirectory="/home/ldelisle/softwares/scriptsForRodriguezCarballoEtAl2020/"
pathForBBCFutils="/home/ldelisle/softwares/" # BBCFutils is found at git clone https://github.com/bbcf/bbcfutils.git -b standalone
mutantGenomeDirectory="/scratch/ldelisle/mutantGenomes/"
FourCDirectory="/scratch/ldelisle/4C_likeHTSstation/"
FourCinRootDirectory="/scratch/ldelisle/all4Cin/"
plotDirectory="/scratch/ldelisle/RC_2020/"
# The pathForChIP should contain the bedgraph from GEO of the ChIP
pathForChIP="/scratch/ldelisle/ChIP/"

# Define the iterations of 4Cin
exclusions="0 1 2 3"
letters="a b c d e"

# Check that the sra table is available:
if [ ! -e ${gitHubDirectory}/tables/sraTable.txt ]; then
  echo "sraTable.txt does not exists. If the paper is published, please send an email to lucille.delisle@epfl.ch"
  exit 1
fi


# First we need to create all the directories
mkdir -p ${mutantGenomeDirectory}
mkdir -p ${FourCDirectory}
for exclusion in ${exclusions}; do
  for letter in "" ${letters}; do
    mkdir -p ${FourCinRootDirectory}4Cin_${exclusion}kb${letter}
  done
done
mkdir -p ${plotDirectory}

# We need seqtk:
if [ `which seqtk | wc -l` -ne 1 ]; then
  echo "When looking for seqtk, got"
  which seqtk
  exit 1
fi

# To prepare the mutant genome creation
# Get the mutant chr2
wget "https://zenodo.org/record/3826913/files/chr2_invTDOM.fa.gz?download=1" -O ${mutantGenomeDirectory}/chr2_invTDOM.fa.gz
gunzip ${mutantGenomeDirectory}/chr2_invTDOM.fa.gz
# Get all chrs from mm10 ucsc
for i in {1..19} X Y M; do
  wget "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr${i}.fa.gz" -O ${mutantGenomeDirectory}/chr${i}.fa.gz
done
# Check the download
# Download the md5sums
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/md5sum.txt -O ${mutantGenomeDirectory}/md5sum.txt
# Get only the one downloaded
cat ${mutantGenomeDirectory}/md5sum.txt | grep -v -P "random|chrUn" > ${mutantGenomeDirectory}/md5sum_mychr.txt
# Change dir to do md5sum
currentDir="$PWD"
cd ${mutantGenomeDirectory}
n=`md5sum -c md5sum_mychr.txt | grep -v OK | wc -l`
cd ${currentDir}

# Check it did not failed
if [ $n -gt 0 ]; then
  echo "Download of chr from mm10 failed"
  exit 1
fi

# Concatenate and make them capital letter length 60
for i in 1 {3..19} X Y M; do
  cat ${mutantGenomeDirectory}/chr${i}.fa.gz >> ${mutantGenomeDirectory}/allChrsExceptchr2_fromUCSC.fa.gz
  echo "chr${i}" >> ${mutantGenomeDirectory}/listOfChrs.txt
done
echo "chr2" > ${mutantGenomeDirectory}/listOfChr2.txt
# This can be long:
seqtk seq -U ${mutantGenomeDirectory}/allChrsExceptchr2_fromUCSC.fa.gz | seqtk subseq -l 60 - ${mutantGenomeDirectory}/listOfChrs.txt | gzip > ${mutantGenomeDirectory}/allChrsExceptchr2.fa.gz
seqtk seq -U ${mutantGenomeDirectory}/chr2.fa.gz | seqtk subseq -l 60 - ${mutantGenomeDirectory}/listOfChr2.txt > ${mutantGenomeDirectory}/chr2_Wt.fa

# You also need to put mm10_rmsk.bed.gz which can be obtained through UCSC website in tools table browser variations and repeat.
if [ ! -e ${mutantGenomeDirectory}/mm10_rmsk.bed.gz ]; then
  echo "You don't have mm10_rmsk.bed.gz you need to download it through the table browser in UCSC"
  exit 1
fi

# To prepare the mapseq step:
# Update the init_bein_minilims.py
cp $gitHubDirectory/scripts/init_bein_minilims.py ${FourCDirectory}
sed -i "s#/mnt/BBCF-raw#${FourCDirectory}#g" ${FourCDirectory}/init_bein_minilims.py
sed -i "s#/home/leleu/htsstation#${pathForBBCFutils}#g" ${FourCDirectory}/init_bein_minilims.py

# To prepare the 4Cin step
for exclusion in ${exclusions}; do
  cp $gitHubDirectory/tables/mm10_viewpointPos_2col.txt ${FourCinRootDirectory}4Cin_${exclusion}kb/
done

# Begin to launch jobs:
# All genomes
jidMakeGenome=$(sbatch --chdir $mutantGenomeDirectory $gitHubDirectory/slurm_scripts/00_makeGenomeFor4CHTS.sh $gitHubDirectory | awk '{print $NF}')

# Get fastq from sra
jidGetFastq=$(sbatch --chdir $FourCDirectory $gitHubDirectory/slurm_scripts/01_getDemultiplexedFastq.sh $gitHubDirectory | awk '{print $NF}')

# Launch the mapping
jidMapping4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidMakeGenome},${jidGetFastq} $gitHubDirectory/slurm_scripts/02_mapping.sh  $gitHubDirectory $mutantGenomeDirectory | awk '{print $NF}')

# Prepare all files for the 4C
jidpre4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidMapping4C} $gitHubDirectory/slurm_scripts/03_prepare4C.sh  $gitHubDirectory $mutantGenomeDirectory | awk '{print $NF}')

# Launch the 4C
jid4C=$(sbatch --chdir $FourCDirectory --dependency afterok:${jidpre4C} $gitHubDirectory/slurm_scripts/04_4Cseq.sh  $gitHubDirectory | awk '{print $NF}')

jid4Cins=""

# Launch all 4Cin
for exclusion in ${exclusions}; do
  jidpre4Cin=$(sbatch --chdir ${FourCinRootDirectory}4Cin_${exclusion}kb --dependency afterok:${jid4C} $gitHubDirectory/slurm_scripts/05_prepare4Cin_smooth.sh  $gitHubDirectory $FourCDirectory $exclusion | awk '{print $NF}')
  for letter in ${letters}; do
    jid4Cin=$(sbatch --chdir ${FourCinRootDirectory}4Cin_${exclusion}kb${letter} --dependency afterok:${jidpre4Cin} $gitHubDirectory/slurm_scripts/06_4Cin.sh  $gitHubDirectory | awk '{print $NF}')
    if [ -z $jid4Cins ]; then
      jid4Cins="${jid4Cin}"
    else
      jid4Cins="${jid4Cins},${jid4Cin}"
    fi
  done
done

# Maybe some 4Cin will fail, you may need to relaunch manually some of them. Do not forget to update the jid4Cins

# Make the averages
jid4CinAv=$(sbatch --chdir ${FourCinRootDirectory} --dependency afterok:${jid4Cins} $gitHubDirectory/slurm_scripts/07_4Cin_average.sh $gitHubDirectory $plotDirectory | awk '{print $NF}')

# Finally do the plots
# The command lines used to build E9_CTCF.bed is in $gitHubDirectory/scripts/CTCF.sh
sbatch --chdir $plotDirectory --dependency afterok:${jid4CinAv} $gitHubDirectory/slurm_scripts/08_plots.sh  $gitHubDirectory $mutantGenomeDirectory $FourCDirectory ${FourCinRootDirectory} $pathForChIP
