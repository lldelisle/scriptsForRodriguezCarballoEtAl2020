mutantGenomeDirectory="/scratch/ldelisle/test/mutantGenomes/"

# The summits of E9_CTCF are extended 100bp each direction
cat E9_CTCF.narrowPeak | awk -v OFS="\t" '$1=="chr2" && $2 < 75800320 && $3 > 73700000 {print $1, $2 + $10 - 100, $2 + $10 + 100, $4}' > E9_CTCF_plotted_extended.bed
# The fasta is extracted
bedtools getfasta -fi ${mutantGenomeDirectory}/chr2_Wt.fa -bed E9_CTCF_plotted_extended.bed -name > E9_CTCF_plotted_extended.fa
# On the website of http://insulatordb.uthsc.edu/ in CTCFBS Prediction Tool the fasta is uploaded.
# The table output for MIT_LM7 is copied into a file:
OUTPUT=E9_CTCF_plotted_extended_insdb_output.txt
# The orientation of the CTCF motif is deduced
cat $OUTPUT | awk -v OFS="\t" '{split($3, a, ":|-"); print a[1], a[2] + $4 + $5/2, a[2] + $4 + $5/2 + 1, $3, $7, $6}' > E9_CTCF.bed
