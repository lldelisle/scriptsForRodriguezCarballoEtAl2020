brFiles=$1
primerFile=$2
pathWithBRFiles=$3
pathForScripts=$4

brFilesSpace=`echo $brFiles | tr "," " "`
for f in $brFilesSpace; do
  ln -s ${pathWithBRFiles}$f .
done

# Extract from the 4cfile the position of the viewpoint
# and put it into a bed
cat $primerFile | awk -F "|" -v OFS="\t" 'NR%2==1{split($3,a,":|-");print a[1],a[2],a[3],$1,1,"+"}' > ${primerFile}_vp.bed

# We shift the bed 
bash $pathForScripts/shiftBedWithMultipleBR.sh $brFiles  ${primerFile}_vp.bed ${primerFile}_vp_shifted.bed 1

mv $primerFile ${primerFile}_ori
# The primer file is formed of groups of 2 lines.
# The first line is description
# The second line is sequence
# We only modify the first lines
h=1
while read line; do
  if [ $h = 1 ]; then
    # We look into the shifted bed which line corresponds to the viewpoint
    # Then we correct the line with the new viewpoint position
    # We also correct the exclude part using the same distance on each side
    cat ${primerFile}_vp_shifted.bed | awk -v l=$line '
      BEGIN{
        split(l,ls,"|")
        written=0
      }
      $4==ls[1]&&written==0{
        oldVP=ls[3]
        ls[3]=$1":"$2"-"$3
        for(i=5;i<=length(ls);i++){
          if(gsub("Exclude=","",ls[i])>0){
            split(ls[i],lse,":|-")
            split(oldVP,lsvp,":|-")
            if($6=="+"){
              ls[i]="Exclude="$1":"$2-(lsvp[2]-lse[2])"-"$3+(lse[3]-lsvp[3])
            }else{
              ls[i]="Exclude="$1":"$2-(lse[3]-lsvp[3])"-"$3+(lsvp[2]-lse[2])
            }
          }
        }
        printf("%s",ls[1])
        for(i=2;i<=length(ls);i++){
          printf("|%s",ls[i])
        }
        print ""
        written=1
      }'  >> $primerFile
    h=0
  else
    echo $line >> $primerFile
    h=1
  fi
done < ${primerFile}_ori
