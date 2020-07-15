##This script takes as input a br file and an annotation file (bed, gtf...)
##It will create where the br file is, one folder per mutant genome with the shifted annotations inside.
##If an annotation overlap an insertion/deletion, it will be split in 2 annotations with the same names.
##Idem if it overlap an inversion.

###############################################
##The script assumes the chromosome is chr2 but you can change it here:
chrToShift<-"chr2"
###############################################

options(scipen=999)
options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
library(tools)

if (length(commandArgs(TRUE))==0){
  script.basename<-getSrcDirectory(function(x) {x})
  cat("Choose the br file.\n")
  pathForBr<-file.choose()
  cat("Choose the file to convert.\n")
  fileToConvert<-file.choose()
  cat("Which is the number of the column with the chromosme name.\n")
  colChr<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the number of the column with the start position.\n")
  colStart<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the number of the column with the end position.\n")
  colEnd<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the number of the column with the strand information (put 0 if there is no).\n")
  colStrand<-as.numeric(readLines(con=stdin(),n=1))
  cat("Which is the name you want to add to each output.\n")
  outputName<-readLines(con=stdin(),n=1)
} else{
  if (commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript shiftAnnot_splitOV.R pathForBr fileToConvert colCr colStart colEnd [colStrand] outputName \n")
    stop()
  }
  pathForBr<-commandArgs(TRUE)[1]
  fileToConvert<-commandArgs(TRUE)[2]
  colChr<-as.numeric(commandArgs(TRUE)[3])
  colStart<-as.numeric(commandArgs(TRUE)[4])
  colEnd<-as.numeric(commandArgs(TRUE)[5])
  if (length(commandArgs(TRUE))>6){
    colStrand<-as.numeric(commandArgs(TRUE)[6])
    outputName<-commandArgs(TRUE)[7]
  } else{
    colStrand<-0
    outputName<-commandArgs(TRUE)[6]
  }
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.basename <- dirname(script.name)
}
################################################################################
# Source the functions which are in the other file:
other.name <- file.path(script.basename, "shiftAnnotFunctions.R")
source(other.name)

cat("Loading input file...")
annotationData <- usefulLDfunctions:::.readFileFromConditionOnNcols(fileToConvert,
                                                                    paste0(">=", max(colChr, colStart, colEnd, colStrand)),
                                                                    keepQuote=T)
cat("loaded.\n")
outputFolder <- dirname(pathForBr)

# Convert to UCSC format
if (!(grepl("chr", annotationData[1, colChr]))){
  annotationData[, colChr] <- paste0(rep("chr", nrow(annotationData)), annotationData[, colChr])
  annotationData[annotationData[, colChr] == "chrMT", colChr] <- "chrM"
  todelete <- grep("^chr(GL|JH)", annotationData[, colChr])
  annotationData <- annotationData[-todelete, ]
}
# Load the br file
brdf <- read.delim(pathForBr, stringsAsFactors=F)
brdf[brdf == ""] <- NA

# Shift the fileToConvert for all the genomes of the br file.
for (i in 1:nrow(brdf)){
  cat("#####################\n")
  print(brdf[i, ])
  cat("####################\n")
  genome <- brdf$genome[i]
  shiftedAnnotations <- shiftDFFromBR(annotationData, genome, brdf, colChr, colStart, colEnd,
                                      colStrand, verbose=T, chromoWithTg=chrToShift, splitIfOverlap=T)
  cat("\n\n\n\n")
  dir.create(file.path(outputFolder, genome), showWarnings=F)
  write.table(shiftedAnnotations, paste0(outputFolder, "/", genome, "/", genome, "_vSplit_", outputName),
              sep='\t', row.names=F, col.names=F, quote=F)
}
