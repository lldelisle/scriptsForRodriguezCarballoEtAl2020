options(scipen=999)
options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
rm(list=ls())

###### Install missing packages and load them ######
safelyLoadAPackageInCRANorBioconductor("zoo")

#### GET PARAMETER FILE ####
if (length(commandArgs(TRUE)) > 0){
  f <- commandArgs(TRUE)[1]
} else{
  # Ask for the config file
  f <- file.choose()
}


#### CHECK EVERYTHING IS READY ####
# Check the parameter file exists and source it
if (! file.exists(f)){
  stop("This file does not exist.")
}
source(f)

# Check if the 4C_function is defined and present
if (! exists("pathForFunctions4C")){
  stop("No path for functions.")
} else if (! file.exists(pathForFunctions4C)){
  stop(paste("This file does not exists:", pathForFunctions4C))
}
source(pathForFunctions4C)

# Check the required parameters
pathForVPNameAndPos <- checkFile("pathForVPNameAndPos")
vpPosDF <- read.delim(pathForVPNameAndPos, h=F)
if (! ncol(vpPosDF) == 2){
  stop("Invalid VPNameAndPos table, expected 2 columns")
}
if (!is.numeric(vpPosDF[1, 2])){
  stop("Invalid VPNameAndPos table, expected numbers only in column 2 (no header)")
}
if( !exists("myChr")){
  stop("myChr is not defined but required to output bedgraphs.")
}
folderWithInputs <- checkDirectory("folderWithInputs")

#### PREPARE INPUT FROM 4C BW ####
# Find all files which looks like segToFrag_BLABLABLA.bw in the folder specified
allBW <- list.files(folderWithInputs, "segToFrag_.*\\.bw$")
# If you have, check if the inputMat are already present.
if (length(allBW) > 0){
  # Guess the name of the corresponding inputMat
  allInputMatOfBW <- paste0(sapply(allBW, function(s){gsub(".bw$|^segToFrag_", "", s)}), ".inputMat")
  names(allInputMatOfBW) <- allBW
} else {
  allInputMatOfBW <- character(0)
}
# Find all inputMat
matrices <- list.files(folderWithInputs, ".inputMat$")
# Check if you have bw file which does not have corresponding inputMat
missingInputMat <- setdiff(allInputMatOfBW, matrices)

if (length(missingInputMat) > 0){
  # You will need a template
  template <- checkFile("template")
  # Select the bw which does not have inputMat
  BWToTransform <- names(allInputMatOfBW[allInputMatOfBW %in% missingInputMat])
  # If rtracklayer is not installed, install it.
  safelyLoadAPackageInCRANorBioconductor("rtracklayer")
  cat("Reading template segToFrag...")
  templateDF <- readBed(template)[, 1:3]
  cat("Done.\n")
  # Restrict it to the chromosome of interest
  templateDF <- subset(templateDF, chr == myChr)
  if(nrow(templateDF) == 0){
    stop("myChr is not part of the template you provided.")
  }
  # put the folder as prefix to be able to generate inputMat at the same place
  filesToWriteInput <- lapply(BWToTransform, function(s){paste0(folderWithInputs, "/", s)})
  cat("Writing missing inputMat from segToFrag bigwigs...")
  # This will write inputMat which asign to each middle of fragment 0 if it was not sequenced and the score in the bw if it is present.
  writeInputMatFromBW(filesToWriteInput, templateDF, myChr)
  cat("Done.\n")
}

#### PREPARE INPUT FROM UMI-4C BDG ####
# Find all bedGraphs files in the folder specified
allBDG <- list.files(folderWithInputs, ".b(e)*(d)*[g|G](raph)*(\\.gz)*$")
# If you have, check if the inputMat are already present.
if (length(allBDG) > 0){
  # Guess the name of the corresponding inputMat
  allInputMatOfBDG <- paste0(sapply(allBDG, tools::file_path_sans_ext), ".inputMat")
  names(allInputMatOfBDG) <- allBDG
} else {
  allInputMatOfBDG <- character(0)
}
# Find all inputMat
matrices <- list.files(folderWithInputs, ".inputMat$")
# Check if you have bw file which does not have corresponding inputMat
missingInputMat <- setdiff(allInputMatOfBDG, matrices)

if (length(missingInputMat) > 0){
  # Select the bdg which does not have inputMat
  BDGToTransform <- names(allInputMatOfBDG[allInputMatOfBDG% in% missingInputMat])
  # put the folder as prefix to be able to generate inputMat at the same place
  filesToWriteInput <- lapply(BDGToTransform, function(s){paste0(folderWithInputs, "/", s)})
  cat("Writing missing inputMat...")
  # This will write inputMat.
  writeInputMatFromBedGraph(filesToWriteInput, NULL, myChr)
  cat("Done.\n")
}

##############
matrices <- list.files(folderWithInputs, ".inputMat$")
if (length(matrices) >= 1){
  ###Check the analysis required parameters####
  # For each matrix we should find the vp position
  colnames(vpPosDF) <- c("name", "pos")
  allVpPos <- NULL
  # For each inputMat file, will look for the viewpoint:
  # From all viewpoint name in the viewPointPosition file, it will use the longest which is included in the file name.
  # For example BlablaHoxd13Blabla, both Hoxd1 and Hoxd13 are included so Hoxd13 (higher number of characters) will be chosen.
  # So if you have an ambiguity, put in the viewPointPosition.txt file as name the full file name.
  for(inputMatFile in matrices){
    curVP <- unique(findVPcoo(inputMatFile, vpPosDF))
    if (! length(curVP) == 1){
      if (length(curVP) == 0){
        stop(paste("For", inputMatFile, ": no viewpoint found in the name. Please change the VPNameAndPos file."))
      } else {
        stop(paste("For", inputMatFile, ":", length(vpPos), " viewpoints found in the name with the coordinates:", paste(curVP, collapse=","), ". Please change the VPNameAndPos file."))
      }
    }
    allVpPos <- c(allVpPos, curVP)
  }
  names(allVpPos) <- matrices
  
  # Check the outputFormat and if HTSStyle is part of them check you have the template.
  outputFormat <- unlist(checkStrings("outputFormat", c("TanayStyle", "ClassicStyle", "HTSStyle"), "TanayStyle"))
  if ("HTSStyle" %in% outputFormat){
    # You will need a template
    if (! exists("templateDF")){
      template <- checkFile("template","NA")
      if (template == "NA"){
        cat("HTSStyle was set as outputFormat but no template is provided.\n")
        if (length(outputFormat) == 1){
          stop("This was the only output. The script stops.")
        } else {
          outputFormat <- setdiff(outputFormat, "HTSStyle")
          cat("There will not be HTSStyle output.\n")
        }
      } else {
        cat("Reading template segToFrag...")
        templateDF <- readBed(template)[, 1:3]
        cat("Done.\n")
        # Restrict it to the chromosome of interest
        templateDF <- subset(templateDF, chr == myChr)
        if (nrow(templateDF) == 0){
          cat("HTSStyle was set as outputFormat but myChr is not part of the template you provided.\n")
          if (length(outputFormat) == 1){
            stop("This was the only output. The script stops.")
          } else {
            outputFormat <- setdiff(outputFormat, "HTSStyle")
            cat("There will not be HTSStyle output.\n")
          }
        }
      }
    }
  }
    
  distanceToRemoveAroundVP <- unlist(checkNumericalValues("distanceToRemoveAroundVP", default = 0))[1]
  if (distanceToRemoveAroundVP < 0){
    cat("distanceToRemoveAroundVP is negative. 0 will be used.\n")
    distanceToRemoveAroundVP <- 0
  }
  smoothThrough <- checkLogicalValue("smoothThrough", F)
  smoothAfter <- checkLogicalValue("smoothAfter", F)
  replaceNAby0 <- checkLogicalValue("replaceNAby0", F)
  if (! exists("outputFolder")){
    outputFolder <- paste0(folderWithInputs, "/", gsub(" ", "_", Sys.time()))
    cat("The bedgraphs will be in :")
    cat(outputFolder)
    cat("\n")
  } else{
    if (! dir.exists(outputFolder)){
      dir.create(outputFolder)
    }
  }
  
  nbsOfFragmentPerWindow <- unlist(checkNumericalValues("nbsOfFragmentPerWindow", 11))
  
  # Add the folder as prefix
  matrices <- sapply(matrices, function(s){paste0(folderWithInputs, "/", s)})
  
  
  ### DO THE SMOOTHING ###
  
  for(nbOfFragmentPerWindow in nbsOfFragmentPerWindow){
    cat("nbOfFragmentPerWindow=", nbOfFragmentPerWindow, "\n")
    for(i in 1:length(matrices)){
      inputMatFile <- matrices[i]
      name <- head(strsplit(basename(tools::file_path_sans_ext(inputMatFile)), "_rep")[[1]], 1)
      vp <- allVpPos[i]
      profile <- read.delim(inputMatFile, h=F)
      colnames(profile) <- c("pos", "score")
      if(anyDuplicated(profile$pos) != 0){
        profile <- aggregate(list(score=profile$score), by=list(pos=profile$pos), FUN=sum)
      }
      profile <- profile[order(profile$pos), ]
      profile$chr <- myChr
      if(smoothAfter){
        profile$score[abs(profile$pos - vp) < distanceToRemoveAroundVP] <- NA
      }
      if(smoothThrough){
        smooth <- profile
        smooth$smoothed <- rollMean(smooth$score, nbOfFragmentPerWindow)
      } else {
        before <- profile[profile$pos < vp, ]
        after <- profile[profile$pos > vp, ]
        before$smoothed <- rollMean(before$score, nbOfFragmentPerWindow)
        after$smoothed <- rollMean(after$score, nbOfFragmentPerWindow)
        smooth <- rbind(before, after)
      }
      if (! smoothAfter){
        smooth$score[abs(smooth$pos - vp) < distanceToRemoveAroundVP] <- NA
      }
      if (replaceNAby0){
        smooth$score[is.na(smooth$score)] <- 0
      }
      if ("TanayStyle" %in% outputFormat){
        smooth$smoothedPos <- round(.p4cWinGeoMeanCoordinate(smooth$pos, nbOfFragmentPerWindow, vp))
        smooth$smoothedPosEnd <- smooth$smoothedPos + 1
        header <- paste0("track type=bedGraph name=", name, "_smoothed", nbOfFragmentPerWindow, "frags_TanayStyle\n")
        outputFile <- paste0(outputFolder, "/", name, "_smoothed", nbOfFragmentPerWindow, "frags_TanayStyle.bedGraph")
        cat(header, file=outputFile)
        temp.df <- smooth[! is.na(smooth$smoothed), c("chr", "smoothedPos", "smoothedPosEnd", "smoothed")]
        temp.df <- temp.df[order(temp.df$smoothedPos), ]
        write.table(temp.df, file=outputFile,
                    sep="\t", row.names = F, col.names = F, quote = F, append = T)
      }
      if ("ClassicStyle" %in% outputFormat){
        smooth$posEnd <- smooth$pos + 1
        header <- paste0("track type=bedGraph name=", name, "_smoothed", nbOfFragmentPerWindow, "frags_classicalStyle\n")
        outputFile <- paste0(outputFolder, "/", name, "_smoothed", nbOfFragmentPerWindow, "frags_classicalStyle.bedGraph")
        cat(header, file=outputFile)
        temp.df <- smooth[!is.na(smooth$smoothed), c("chr", "pos", "posEnd", "smoothed")]
        temp.df <- temp.df[order(temp.df$pos), ]
        write.table(temp.df, file=outputFile,
                    sep="\t", row.names = F, col.names = F, quote = F, append = T)
      }
      if("HTSStyle" %in% outputFormat){
        header <- paste0("track type=bedGraph name=", name, "_smoothed", nbOfFragmentPerWindow, "frags_HTSStyle\n")
        outputFile <- paste0(outputFolder, "/", name, "_smoothed", nbOfFragmentPerWindow, "frags_HTSStyle.bedGraph")
        cat(header, file=outputFile)
        templateDF$pos <- round((templateDF$start + templateDF$end) / 2)
        mdf <- merge(templateDF, smooth, all=T)
        temp.df <- mdf[!is.na(mdf$smoothed) & !is.na(mdf$start), c("chr", "start", "end", "smoothed")]
        temp.df <- temp.df[order(temp.df$start), ]
        write.table(temp.df, file=outputFile,
                    sep="\t", row.names = F, col.names = F, quote = F, append = T)
      }
    }
  }
} else {
  cat("There are no matrice to smooth.\n")
}
