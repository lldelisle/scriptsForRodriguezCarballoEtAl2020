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
pathWithResultsOf4C=$3
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

exists=`conda info --envs | awk '$1=="pgt_3.5"{print}' | wc -l`
if [ $exists -ne 1 ]; then
  conda env create -n pgt_3.5 pygenometracks=3.5  
fi
conda activate pgt_3.5

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
bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} ${gitHubDirectory}/tables/subTADs.bed subTADs_oninvTDOM.bed 0
sort -k1,1 -k2,2n subTADs_oninvTDOM.bed > subTADs_oninvTDOM_sorted.bed
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
bash ${pathForScripts}/generateAllIni.sh ${pathWithIniFiles} ${gitHubDirectory} ${FourCinRootDirectory} ${pathWithResultsOf4C} ${pathForChIP}

# Plot

# Figure 1A
pgt --tracks ${pathWithIniFiles}/fig1A.ini --region chr2:73700000-75800000 -o fig1A.pdf

# Figure 1B
pgt --tracks ${pathWithIniFiles}/fig1B.ini --region chr2:74400001-75800000 -o fig1B.pdf

for i in 2 3 4 5A 5C S2B S3AB S4BC S5B S6AB S6C S7A S7B; do
  pgt --tracks ${pathWithIniFiles}/fig${i}.ini --region chr2:74401941-75800320 -o fig${i}.pdf
done

# Figure S2A
# UCSC chr2:75131700-75156173

# Figure S3C
# UCSC

# Figure S5A
# UCSC chr2:75131700-75156173

# For the quantification:
# Peaks need to be called from E9_wt_Hoxd11, Hoxd9 and Hoxd4
mkdir -p callPeaks
for hoxd in 4 9 11; do
  ln -s ${pathWithResultsOf4C}/toGEO/segToFrag_E9_FLB_wt_Hoxd${hoxd}.bw callPeaks/
done

# Prepare a config file for peak calling
template=mm10_NlaIIIDpnII.bedGraph.gz
echo "###Required parameters
coordinatesToPlot <- c(\"chr2:74401941-75800320\")
pathForFunctions4C <- \"${pathForScripts}/4C_functions.R\"
pathForVPNameAndPos <- \"$gitHubDirectory/tables/mm10_viewpointPos_2col.txt\"
folderWithInputs <- \"$PWD/callPeaks/\"

###Optional parameters
#Output location
outputPath <- \"$PWD/callPeaks/plots\"
#Output format
usePng <- F #If you want to use png replace F by T
pngRes <- 96 #This is the resolution of the png file.
#Additional features on the plot
plotModel <- T #Do you want to see a blue line corresponding to the model to which the score will be compared.

##Analysis parameters:
wins <- 1e6 #Relatively to the viewpoint to which extend do you want to look for peaks
qWrValues <- 0.6 # Parameter on the ratio over background
qWdValues <- 1.5 # Parameter on the difference with background
nbsOfFragmentPerWindow <- 11 # smoothing window

### 4C PARAMETERS ###
template <- \"${template}\"
" > configFile_for_peak_calling.R

# Generate the template
possibleT=`ls ${pathWithResultsOf4C}/res_files_4Cseq_Wt*/normalised_scorePerFeature_*.bedGraph.gz | awk 'NR==1{print $1}'`
if [ ! -z ${possibleT} ]; then
  zcat ${possibleT} | awk 'NF==4{print $1"\t"$2"\t"$3}' | gzip > ${template}
else
  echo "There is no file ${genome}/normalised_scorePerFeature_*.bedGraph.gz"
  exit 1
fi

if [ ! -e ${pathForScripts}/callPeaks.R ]; then
  echo "${pathForScripts}/callPeaks.R does not exists"
  exit 1
fi

Rscript ${pathForScripts}/callPeaks.R configFile_for_peak_calling.R

# First the CS93, CS65, CS39
cat ${gitHubDirectory}/tables/annotations_enhancers_figure.bed | grep -P "CS93|CS65|CS39" > enhancers.bed
peakCalling_file=callPeaks/E9_FLB_wt_Hoxd11.inputMat_74685562_1_11_1.5_0.6_peaks.bed
bedtools intersect -a $peakCalling_file -b enhancers.bed -wa -wb | cut -f 1-3,7 > quantification.bed

# Then ELCR2
cat ${gitHubDirectory}/tables/annotations_enhancers_figure.bed | grep -w "ELCR2" > enhancers.bed
peakCalling_file=callPeaks/E9_FLB_wt_Hoxd9.inputMat_74697247_1_11_1.5_0.6_peaks.bed
bedtools intersect -a $peakCalling_file -b enhancers.bed -wa -wb | cut -f 1-3,7 >> quantification.bed

# Then ELCR3
cat ${gitHubDirectory}/tables/annotations_enhancers_figure.bed | grep -w "ELCR3" > enhancers.bed
peakCalling_file=callPeaks/E9_FLB_wt_Hoxd4.inputMat_74723268_1_11_1.5_0.6_peaks.bed
bedtools intersect -a $peakCalling_file -b enhancers.bed -wa -wb | cut -f 1-3,7 >> quantification.bed

# For HoxD: From the end of Prox to the TSS of Hoxd1:
left_border=`cat ${gitHubDirectory}/tables/annotations_enhancers_figure.bed | grep Prox | cut -f 3`
cat genes_around_TDOM_customized.bed | awk -v OFS="\t" -v lb=$left_border '$4=="Hoxd1"{$2=lb;print $1,$2,$3,"HoxD"}' >> quantification.bed

# For TDOM we use the subTADs
cat ${gitHubDirectory}/tables/subTADs.bed | awk -v OFS="\t" 'NR==2{start=$2}END{$2=start;print $1,$2,$3,"TDOM"}' >> quantification.bed

# We need to add the same for invTDOM:
bash ${pathForScripts}/shiftBedWithMultipleBR.sh ${brFiles} quantification.bed quantification_oninvTDOM.bed 0
cat quantification_oninvTDOM.bed | awk -v OFS="\t" '{$4=$4"_oninvTDOM";print}' > all_quantification.bed
cat quantification.bed >> all_quantification.bed

# We launch the python script which will compute statistics on 4C
# Quantification on virtual 4C
# Quantification on 4C
python ${gitHubDirectory}/scripts/stats.py ${gitHubDirectory}/tables ${FourCinRootDirectory} ${pathWithResultsOf4C}/toGEO/ \
  all_quantification.bed ${gitHubDirectory}/tables/quantifications_4C_table.txt stats_4c.txt quantif_4c.txt > quantif_vHiC.txt
