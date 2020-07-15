# Generate a h5 corrected matrix from a validPair from Dekker pipeline
# Require HiCExplorer
validPairs=$1
sizeFile=$2
binSize=$3
outputFile=$4
# We create a temporary directory
TEMPORARYDIR=`mktemp -d -t`
if [[ ${validPairs} =~ \.gz$ ]]; then
  gunzip -c ${validPairs} > ${TEMPORARYDIR}/input.txt
else
  if [[ ! "${validPairs}" = "/"* ]]; then
    validPairs="$PWD/${validPairs}"
  fi
  ln -s ${validPairs} ${TEMPORARYDIR}/input.txt
fi
# We remove duplicates
awk '{
  if(chrB1 != $2 || chrB2 != $8 || posB1 != $3 || posB2 != $9){
    print
  }
  chrB1=$2
  chrB2=$8
  posB1=$3
  posB2=$9
}' ${TEMPORARYDIR}/input.txt > ${TEMPORARYDIR}/unique.txt

cooler csort -i tabix -c1 2 -c2 8 -p1 3 -p2 9 -o ${TEMPORARYDIR}/unique.csort.gz ${TEMPORARYDIR}/unique.txt ${sizeFile}

cooler makebins -o ${TEMPORARYDIR}/bins.txt ${sizeFile} ${binSize}

cooler cload tabix --assembly mm10 -c2 8 -p2 9 ${TEMPORARYDIR}/bins.txt ${TEMPORARYDIR}/unique.csort.gz ${TEMPORARYDIR}/raw.cool

hicCorrectMatrix correct --matrix ${TEMPORARYDIR}/raw.cool --iterNum 500 --outFileName ${outputFile} --filterThreshold -1.5 5.0
