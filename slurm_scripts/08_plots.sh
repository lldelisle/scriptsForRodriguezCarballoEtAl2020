#!/bin/bash -l

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 50G
#SBATCH --time 6:00:00
#SBATCH --job-name plots

gitHubDirectory=$1
mutantGenomeDirectory=$2
FourCDirectory=$3
FourCinRootDirectory=$4
pathForChIP=$5

module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0

path="$PWD/"
pathForScripts="${gitHubDirectory}/scripts/"
pathWithTableWithGenomes="${gitHubDirectory}/tables/table.txt"
pathWithBRFiles="${gitHubDirectory}/tables/"
pathWithIniFiles="${gitHubDirectory}/ini_files/"
pathForSizes="${mutantGenomeDirectory}/Wt/mm10.fa.fai"
pathWithIniFiles="${path}/"

# For the first figure:

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

# The Hi-C matrix was obtained with an old version of hicexplorer:
conda create -n hicexplorer2.1.4 --yes hicexplorer=2.1.4 python=2.7 tabix
conda activate hicexplorer2.1.4

# Create the matrix file for figure 1
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2713717&format=file&file=GSM2713717%5FPL%5FHiC%5FE12%5FWt%2EvalidPair%2Etxt%2Egz" -O E12_PFL_HiC_validPairs.txt.gz

bash ${pathForScripts}/generateH5fileFromValidPairs.sh E12_PFL_HiC_validPairs.txt.gz ${pathForSizes} 40000 E12_PFL_HiC.h5

conda deactivate

# Download the CTCF coverage for Fig1A
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM2713708&format=file&file=GSM2713708%5FPFL%5FE12%5FWt%5FCTCF%2EbedGraph%2Egz" -O PFL_E12_Wt_CTCF.bedGraph.gz
gunzip PFL_E12_Wt_CTCF.bedGraph.gz
grep "chr2" PFL_E12_Wt_CTCF.bedGraph > PFL_E12_Wt_CTCF_chr2.bedGraph
bedGraphToBigWig PFL_E12_Wt_CTCF_chr2.bedGraph ${pathForSizes} PFL_E12_Wt_CTCF_chr2.bw

# Make the bigwig for ChIP if you got the bedgraph from GEO:
for f in ${pathForChIP}/*.bedGraph.gz; do
  name=`basename $f .bedGraph.gz | awk '{split($1, a, "_E9"); if(length(a) == 2){print "E9"a[2]}else{print $0}}'`
  if [ ! -e ${pathForChIP}/${name}.bw ]; then
    gunzip -c $f | grep chr2 > temp.bdg
    bedGraphToBigWig temp.bdg ${pathForSizes} ${pathForChIP}/${name}.bw
    rm temp.bdg
  fi
done


# For all figures:

# Create a bed with rgb field corresponding to motif orientation of CTCF:
awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" '{if($6 == "+"){color = colorPos}else{color = colorNeg};print $1, $2, $3, $4, $5, $6, $2, $2, color}' ${gitHubDirectory}/tables/E9_CTCF.bed > CTCF_colored.bed

# Create a small gtf for genes:
wget "https://zenodo.org/record/3820860/files/mergeGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.93_ExonsOnly_UCSC.gtf.gz?download=1" -O all_genes.gtf.gz

zcat all_genes.gtf.gz | awk '$1=="chr2" && $4 < 77000000 && $5 > 73000000{print}' > genes_around_TDOM.gtf

cat genes_around_TDOM.gtf | awk -F "\t" '$9~/protein_coding/ || $9~/Haglr/ || $9~/Hoxd3os/{if(!($9~/Hoxd3-203/ || $9~/Hoxd3-os1-201/ || $9~/Haglr-202/)){print}}' > genes_around_TDOM_customized.gtf

cat genes_around_TDOM.gtf | awk -F "\t" '$9~/Gm13652/ || $9~/2600014E21Rik/ {print}' > Hog_Tog.gtf

exists=`conda info --envs | awk '$1=="pgt_rc2020"{print}' | wc -l`
if [ $exists -ne 1 ]; then
  conda env create -f ${pathForScripts}/environment.yml
fi
conda activate pgt_rc2020

# Convert gtf to bed
python ${pathForScripts}/fromgtfTobed12.py --output genes_around_TDOM.bed --mergeTranscriptsAndOverlappingExons genes_around_TDOM.gtf

python ${pathForScripts}/fromgtfTobed12.py --output genes_around_TDOM_customized.bed --mergeTranscriptsAndOverlappingExons genes_around_TDOM_customized.gtf

python ${pathForScripts}/fromgtfTobed12.py --output Hog_Tog.bed --mergeTranscriptsAndOverlappingExons Hog_Tog.gtf

cat genes_around_TDOM_customized.bed | awk -v OFS="\t" '$6=="+"{print $1, $2, $3, $4}' > genes_around_TDOM_customized_pos.bed

cat genes_around_TDOM_customized.bed | awk -v OFS="\t" '$6=="-"{print $1, $2, $3, $4}' > genes_around_TDOM_customized_neg.bed

# Shift all annotations to invTDOM
# Get the brFiles from the table
if [ ! -e ${pathWithTableWithGenomes} ]; then
  echo "${pathWithTableWithGenomes} does not exists"
  exit 1
fi

brFiles=`cat ${pathWithTableWithGenomes} | awk '$1=="invTDOM"{print $2}'`

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
bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} ${gitHubDirectory}/tables/E9_CTCF.bed CTCF_oninvTDOM.bed 1
bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} ${gitHubDirectory}/tables/Bd.bed Bd_oninvTDOM.bed 0
bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} ${gitHubDirectory}/tables/annotations_enhancers_figure.bed annotations_enhancers_figure_oninvTDOM.bed 0
sort -k1,1 -k2,2n annotations_enhancers_figure_oninvTDOM.bed > annotations_enhancers_figure_oninvTDOM_sorted.bed
if [ ! -e ${pathForScripts}/shiftGtfWithMultipleBR.sh ]; then
  echo "${pathForScripts}/shiftGtfWithMultipleBR.sh does not exists"
  exit 1
fi
bash ${pathForScripts}/shiftGtfWithMultipleBR.sh ${brFiles} ${path}genes_around_TDOM_customized.gtf ${path}genes_around_TDOM_customized_oninvTDOM.gtf 
bash ${pathForScripts}/shiftGtfWithMultipleBR.sh ${brFiles} ${path}Hog_Tog.gtf ${path}Hog_Tog_oninvTDOM.gtf 

# Convert gtf to bed
python ${pathForScripts}/fromgtfTobed12.py --output ${path}genes_around_TDOM_customized_oninvTDOM.bed --mergeTranscriptsAndOverlappingExons ${path}genes_around_TDOM_customized_oninvTDOM.gtf

python ${pathForScripts}/fromgtfTobed12.py --output Hog_Tog_oninvTDOM.bed --mergeTranscriptsAndOverlappingExons Hog_Tog_oninvTDOM.gtf

# Create a bed with rgb field corresponding to motif orientation of CTCF:
sort -k1,1 -k2,2n CTCF_oninvTDOM.bed | awk -F "\t" -v OFS="\t" -v colorPos="236,28,36" -v colorNeg="46,49,145" '{if($6 == "+"){color = colorPos}else{color = colorNeg};print $1, $2, $3, $4, $5, $6, $2, $2, color}' > CTCF_colored_oninvTDOM.bed

cat genes_around_TDOM_customized_oninvTDOM.bed | awk -v OFS="\t" '$6=="+"{print $1, $2, $3, $4}' > genes_around_TDOM_customized_pos_oninvTDOM.bed

cat genes_around_TDOM_customized_oninvTDOM.bed | awk -v OFS="\t" '$6=="-"{print $1, $2, $3, $4}' > genes_around_TDOM_customized_neg_oninvTDOM.bed


# Create the bed file for viewpoints used in 4Cin
while read line; do
  experiment=`echo $line | awk '{print $2}'`
  primerFile=`ls ${FourCinRootDirectory}/*/model_${experiment}_onTDOM/primers.txt | head -n 1`
  bedFile=${experiment}_viewpoints.bed
  while read l; do
    vp=`echo $l | awk '{print $1}'`
    chrom=`echo $l | awk -F " |:" '{print $2}'`
    pos=`echo $l | awk -F " |:" '{print $3}'`
    echo -e "${chrom}\t${pos}\t$((${pos}+1))\t${vp}\t0\t-" >> ${bedFile}
    echo -e "${chrom}\t${pos}\t$((${pos}+1))\t${vp}\t0\t+" >> ${bedFile}
  done < ${primerFile}
  sort -k1,1 -k2,2n ${bedFile} > temp.bed
  echo "track type=bed name=${experiment}_viewpoints" > ${bedFile}
  cat temp.bed >> ${bedFile}
  rm temp.bed
done < <(sed 1d ${gitHubDirectory}/tables/4CinTable.txt)


# Generate all ini files using the paths
bash ${pathForScripts}/generateAllIni.sh ${pathWithIniFiles} ${gitHubDirectory} ${FourCinRootDirectory} ${FourCDirectory} ${pathForChIP}

# Plot

# Figure 1A
pgt --tracks ${pathWithIniFiles}/fig1A.ini --region chr2:73700000-75800000 -o fig1A.pdf

# Figure 1B
pgt --tracks ${pathWithIniFiles}/fig1B.ini --region chr2:74400001-75800000 -o fig1B.pdf

# Figure 2
pgt --tracks ${pathWithIniFiles}/fig2.ini --region chr2:74401941-75800320 -o fig2.pdf

# Figure 3
pgt --tracks ${pathWithIniFiles}/fig3.ini --region chr2:74401941-75800320 -o fig3.pdf

# Figure 4
pgt --tracks ${pathWithIniFiles}/fig4.ini --region chr2:74401941-75800320 -o fig4.pdf

# Figure 5A
pgt --tracks ${pathWithIniFiles}/fig5A.ini --region chr2:74401941-75800320 -o fig5A.pdf

# Figure 5C
pgt --tracks ${pathWithIniFiles}/fig5C.ini --region chr2:74401941-75800320 -o fig5C.pdf

# Figure S2A
# UCSC chr2-75131700-75156173


# Figure S2BC
pgt --tracks ${pathWithIniFiles}/figS2BC.ini --region chr2:74401941-75800320 -o figS2BC.pdf

# Figure S2D
# UCSC chr2:74401941-75800320

# Figure S3BC
pgt --tracks ${pathWithIniFiles}/figS3BC.ini --region chr2:74401941-75800320 -o figS3BC.pdf

# Figure S4A
# UCSC chr2-75131700-75156173


# Figure S5
pgt --tracks ${pathWithIniFiles}/figS5.ini --region chr2:74401941-75800320 -o figS5.pdf

# Figure S6
pgt --tracks ${pathWithIniFiles}/figS6.ini --region chr2:74401941-75800320 -o figS6.pdf
