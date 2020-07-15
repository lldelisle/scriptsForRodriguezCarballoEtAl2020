scriptFullPath=$0
multipleBR=$1
originalFile=$2
finalFile=$3
useStrand=$4

scriptFolder=`dirname $scriptFullPath`
brFiles=`echo $multipleBR | tr "," " "`
fileToShift=$originalFile
for br in $brFiles; do
  if [ $useStrand = 1 ]; then
    Rscript $scriptFolder/shiftAnnot_splitIfOV.R $br $fileToShift 1 2 3 6 .bed
  else
    Rscript $scriptFolder/shiftAnnot_splitIfOV.R $br $fileToShift 1 2 3 .bed
  fi
  folderBR=`dirname $br`
  fileToShift=`cat $br | awk -v f=$folderBR 'NR==2{print f"/"$1"/"$1"_vSplit_.bed"}'`
done
mv $fileToShift $finalFile
  
