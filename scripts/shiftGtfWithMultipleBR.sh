scriptFullPath=$0
multipleBR=$1
originalFile=$2
finalFile=$3

scriptFolder=`dirname $scriptFullPath`
brFiles=`echo $multipleBR | tr "," " "`
fileToShift=$originalFile
for br in $brFiles; do
  Rscript $scriptFolder/shiftAnnot_splitIfOV.R $br $fileToShift 1 4 5 7 .gtf
  folderBR=`dirname $br`
  fileToShift=`cat $br | awk -v f=$folderBR 'NR==2{print f"/"$1"/"$1"_vSplit_.gtf"}'`
done
mv $fileToShift $finalFile
  
