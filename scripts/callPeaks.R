options(scipen=999)
options(stringsAsFactors = F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
rm(list=ls())

###### Install missing packages and load them ######
safelyLoadAPackageInCRANorBioconductor("isotone")
devtools::install_github("deWitLab/peakC")
library(peakC)
safelyLoadAPackageInCRANorBioconductor("caTools")

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

# CoordinatesToPlot
if(!exists("coordinatesToPlot")){
  stop("No coordinatesToPlot defined.")
} 
coordinatesToPlotDF <- as.data.frame(
                         matrix(
                           unlist(
                             lapply(coordinatesToPlot,
                                    function(s){
                                      v <- unlist(
                                             strsplit(as.character(s), ":|-")
                                           )
                                      if (length(v) == 3){
                                        return(v)
                                      } else {
                                        stop("Invalid coordinate, expected chrN:start-end")
                                      }
                                    })
                           ), ncol = 3, byrow = T)
                         )
colnames(coordinatesToPlotDF) <- c("chr", "start", "end")
coordinatesToPlotDF[, c(2, 3)] <- apply(coordinatesToPlotDF[, c(2, 3)], 2, gsub, pattern=",", replacement="") # To be compatible with UCSC
coordinatesToPlotDF[, c(2, 3)] <- tryCatch(apply(coordinatesToPlotDF[, c(2, 3)],
                                                 2, as.numeric),
                                           warning=function(w){stop("invalid coordinates to plot, expected chrN:start-end")})
myChr <- coordinatesToPlotDF$chr[1]
if (! all(c(coordinatesToPlotDF$chr) == myChr)){
  stop("For the moment it is not possible to plot 2 regions with different chromosomes.")
}

## Check the required parameters
pathForVPNameAndPos <- checkFile("pathForVPNameAndPos")
vpPosDF <- read.delim(pathForVPNameAndPos, h=F)
if (! ncol(vpPosDF) == 2){
  stop("Invalid VPNameAndPos table, expected 2 columns")
}
if (!is.numeric(vpPosDF[1, 2])){
  stop("Invalid VPNameAndPos table, expected numbers only in column 2 (no header)")
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
matrices<-list.files(folderWithInputs,".inputMat$")
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

  #restrict the annotation file
  if (exists("bedFileWithAnnotations")){
    if (file.exists(bedFileWithAnnotations)){
      annot <- readBed(bedFileWithAnnotations)
      #We restrict annot to the plotted region
      annot <- annot[annot[, 1] == myChr & annot[, 2] < max(coordinatesToPlotDF[, 3]) & annot[, 3] > min(coordinatesToPlotDF[, 2]), ]
    } else {
      annot <- data.frame()
    }
  } else {
    annot <- data.frame()
  }
  
  if (! exists("outputPath")){
    outputPath <- paste0(folderWithInputs, "/", gsub(" ", "_", Sys.time()))
    cat("The plot will be in :")
    cat(outputPath)
    cat(".pdf\n")
  } else {
    if (! dir.exists(dirname(outputPath))){
      dir.create(dirname(outputPath))
    }
  }
  
  if (exists("usePng")){
    if (! is.logical(usePng)){
      cat("usePng is not boolean and is not taken into account.\n The output will be pdf.\n")
      usePng <- F
    }
  } else {
    usePng <- F
  }
  
  if (usePng){
    if (! exists("pngRes")){
      cat("You need to specify a resolution for png files.(pngRes<-96 or whatever).\n By default it will be 96")
      pngRes <- 96
    } else {
      if (! is.numeric(pngRes)){
        cat("The resolution for png files need to be a number, for example pngRes<-96 or whatever.\n By default it will be 96")
        pngRes <- 96
      }
    }
  }
  
  if (exists("plotModel")){
    if (! is.logical(plotModel)){
      cat("plotModel is not boolean and is not taken into account.\n The model will not be plotted.\n")
      plotModel <- F
    }
  } else{
    plotModel <- F
  }
  
  wins <- unlist(checkNumericalValues("wins", 1e6))
  qWrValues <- unlist(checkNumericalValues("qWrValues", 1))
  qWdValues <- unlist(checkNumericalValues("qWdValues", 1.5))
  nbsOfFragmentPerWindow <- unlist(checkNumericalValues("nbsOfFragmentPerWindow", 11))
  
  
  # Prepare a dataframe where the thresholds will be stored.
  thresholds <- data.frame(file=matrices)
  # Add the folder as prefix
  matrices <- sapply(matrices, function(s){paste0(folderWithInputs, "/", s)})
  
  
  ### DO THE PEAK CALLING ###

  if (! usePng){
    pdf(paste0(outputPath, ".pdf"), height = (length(matrices) + 1) * 2,
        title = basename(f))
  }
  par(mfrow=c(length(matrices) + 1, 1))
  for (win in wins){
    cat("win=", win, "\n")
    for (nbOfFragmentPerWindow in nbsOfFragmentPerWindow){
      cat("nbOfFragmentPerWindow=", nbOfFragmentPerWindow, "\n")
      for (qWdValue in qWdValues){
        cat("qWdValue=", qWdValue, "\n")
        thresholds[, paste0("w=", win / 1e6, "s=", nbOfFragmentPerWindow, "qWd=", qWdValue)] <- NA
        for(qWrValue in qWrValues){
          cat("qWrValue=", qWrValue, "\n")
          thresholds[, paste0("w=", win / 1e6, "s=", nbOfFragmentPerWindow, "qWr=", qWrValue)] <-NA
          allRes <- list()
          for (i in 1:length(matrices)){
            inputMatFile <- matrices[i]
            name <- head(strsplit(basename(tools::file_path_sans_ext(inputMatFile)), "_rep")[[1]], 1)
            vpPos <- allVpPos[i]
            # Read the matrix and perform quality
            testMat <- readMatrix(inputMatFile, win, vp.pos = vpPos)
            # Call peak
            res <- single.analysis(testMat$data, vp.pos = vpPos, wSize=nbOfFragmentPerWindow, qWd = qWdValue, qWr = qWrValue)
            allRes[[i]] <- res
            # Write a bed file with the peaks (concatenate consecutive fragments)
            cat("track name=", name, "_peakCalling_", paste(win / 1e6, nbOfFragmentPerWindow, qWdValue, qWrValue, sep = "_"), "\n", sep = '', file = paste(inputMatFile, vpPos, win / 1e6, nbOfFragmentPerWindow, qWdValue, qWrValue, "peaks.bed", sep = "_"))
            write.table(bedFromPeakAndAllcoo(res$peak,res$dbR[, 1], myChr), paste(inputMatFile, vpPos, win / 1e6, nbOfFragmentPerWindow, qWdValue, qWrValue, "peaks.bed", sep = "_"),
                        row.names = F, col.names = F, quote = F, sep = "\t", append = T)
            # Write the quality statistics
            if (! file.exists(paste0(inputMatFile, "_qual.txt"))){
              cat(names(testMat$quality), sep="\t", file=paste0(inputMatFile, "_qual.txt"))
              cat("\n", file=paste0(inputMatFile, "_qual.txt"), append = T)
              cat(unlist(testMat$quality), sep="\t", file=paste0(inputMatFile, "_qual.txt"), append = T)
            }
            # Put the values of the threshold used in the thresholds dataframe
            thresholds[i, paste0("w=", win / 1e6, "s=", nbOfFragmentPerWindow, "qWr=", qWrValue)] <- peakC:::getThreshold(res$ratios[, 2], qWrValue)
            thresholds[i, paste0("w=", win / 1e6, "s=", nbOfFragmentPerWindow, "qWd=", qWdValue)] <- peakC:::getThreshold(res$deltas[, 2], qWdValue)
          }
          # Do the plot
          if (usePng){
            png(paste0(outputPath, "_w=", win / 1e6, "_s=", nbOfFragmentPerWindow, "_qWr=", qWrValue, "_qWd=", qWdValue, "_",
                       coordinatesToPlotDF$start[k] / 1e6, "-", coordinatesToPlotDF$end[k] / 1e6, ".png"),
                width = 7,height = (length(matrices) + 1) * 2, units="in", res=pngRes)
          }
          for (k in 1:nrow(coordinatesToPlotDF)){
            coo <- c(coordinatesToPlotDF$start[k], coordinatesToPlotDF$end[k])
            par(mar=c(0,2,0,2))
            for (i in 1:length(matrices)){
              inputMatFile <- matrices[i]
              name <- head(strsplit(basename(tools::file_path_sans_ext(inputMatFile)), "_rep")[[1]], 1)
              res <- allRes[[i]]
              # Plot the smoothed 4C-signal
              plot(res$dbR[, 1], res$dbR[, 2], type='h', ylab="4C signal", xlim=coo, xaxt='n')
              if (plotModel){
                # Plot the background model
                lines(res$dbR[, 1], res$dbR[, 3], lwd=2, col='blue')
              }
              # Plot the peaks in red
              points(res$dbR[res$dbR[, 1] %in% res$peak, 1], res$dbR[res$dbR[, 1] %in% res$peak, 2], type='h', col="red")
              # Add the name of the inputMat
              text(x = coo[1], y= max(res$dbR[, 3]) / 2, labels = name, adj=c(0,0))
              if (i == 1){
                # The first time, write the parameters
                text(x=mean(coo), y=max(res$dbR[,3]) * 0.9, labels=paste("qWr =", qWrValue, "qWd =", qWdValue, "smooth=", nbOfFragmentPerWindow, "win=", win))
              }
            }
            # In the last track, put annotations
            par(mar=c(3,2,3,2))
            plot(res$dbR[, 1], res$dbR[, 2], type='n', xlab="chromosomal position", xlim=coo, yaxt='n', ylim=c(0, 1))
            if (nrow(annot) > 0){
              for (k in 1:nrow(annot)){
                start <- annot[k, 2]
                end <- annot[k, 3]
                ypos <- c(0, 0.7, 0.1, 0.8, 0.2, 0.9)[k %% 6 + 1]
                rect(xleft = start, ybottom = 0.5, xright = end, ytop = 0.6, col = "black", border = NA)
                # The text is left aligned to the start of the annotation
                text(x = start, y=ypos, labels = annot[k, 4], adj = c(0, 0))
              }
            }
            if (usePng){
              dev.off()
            }
          }
        }
      }
    }
  }
  if (! usePng){
    dev.off()
  }
  # The thresholds are written next to the pdf file
  write.table(thresholds, file=paste0(outputPath, "_thresholds.txt"), row.names = F, col.names = T, quote = F, sep = "\t")
} else {
  cat("There are no matrices to plot.\n")
}
