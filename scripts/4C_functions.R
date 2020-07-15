writeInputMatFromDFbdgFormat <- function(df, fileOut){
  df$mid <- round((df[, 2] + df[, 3]) / 2)
  write.table(df[, c("mid", colnames(df)[4])], file=fileOut, row.names=F, col.names=F, quote=F, sep="\t")
}

writeInputMatFromBedGraph <- function(listOfFiles, templateDF, chrToRestrict){
  for (file in listOfFiles){
    segToFrag <- readBedGraph(file)
    segToFrag <- segToFrag[segToFrag[,1] %in% chrToRestrict, ]
    if (is.null(templateDF)){
      zVal <- sum(segToFrag[, 4] == 0)
      if (zVal == 0){
        cat("The file", file, "have no 0 values and no template are provided. This might be an error.\n")
      }
      mergeSTF <- segToFrag
    } else {
      mergeSTF <- add0ToSegToFrag(segToFrag, templateDF)
    }
    fileOut <- paste0(gsub("^segToFrag_", "", tools::file_path_sans_ext(file), ignore.case = T), ".inputMat")
    fileOut <- gsub("/segToFrag_", "", fileOut)
    writeInputMatFromDFbdgFormat(mergeSTF, fileOut)
  }
}

bedGraphFromGR <- function(gr){
  bdg <- data.frame(gr)
  bdg$start <- bdg$start - 1
  bdg <- bdg[, c("seqnames", "start", "end", "score")]
  colnames(bdg)[1] <- "chr"
  return(bdg)
}

add0ToSegToFrag <- function(segToFragDF, templateDF){
  colnames(segToFragDF) <- c(colnames(templateDF),"score")
  mergeSTF <- merge(segToFragDF, templateDF, all=T)
  if (nrow(mergeSTF) > (nrow(templateDF) + 10)){
    stop("Invalid template.")
  }
  mergeSTF <- mergeSTF[order(mergeSTF[, 2]), ]
  mergeSTF $ score[is.na(mergeSTF$score)] <- 0
  return(mergeSTF)
}

writeInputMatFromBW <- function(listOfFiles, templateDF, chrToRestrict){
  for (file in listOfFiles){
    segToFragGR <- import.bw(file)
    segToFragGR <- subset(segToFragGR, seqnames %in% chrToRestrict)
    segToFrag <- bedGraphFromGR(segToFragGR)
    mergeSTF <- add0ToSegToFrag(segToFrag, templateDF)
    fileOut <- paste0(gsub(".bw$|^segToFrag_", "", gsub("/segToFrag_", "/", file, ignore.case=T), ignore.case=T), ".inputMat")
    writeInputMatFromDFbdgFormat(mergeSTF, fileOut)
  }
}

findVPcoo <- function(fn, vpPosDF){
  pos1 <- which(sapply(vpPosDF$name, function(s){grepl(s, fn, ignore.case=T)}))
  if (length(pos1) > 1){
    pos1 <- pos1[which.max(sapply(vpPosDF$name[pos1], nchar))]
  }
  return(vpPosDF[pos1, "pos"])
}

vpPosFromBedWithFragmentCoo <- function(bedDF){
  vpPosDF <- data.frame(pos=bedDF[, 2])
  vpPosDF$pos[bedDF[, 6] == "+"] <- bedDF[bedDF[, 6] == "+", 3]
  vpPosDF$name <- bedDF[, 4]
  return(vpPosDF)
}

bedFromPeakAndAllcoo <- function(peakCoo, allCoo, myChr){
  df <- data.frame(index=which(allCoo %in% peakCoo))
  df$pos <- allCoo[df$index]
  df$nextIndex <- c(df$index[-1], NA)
  df$isConsec <- df$nextIndex == (df$index + 1)
  df$isStart <- c(T, !df$isConsec[-nrow(df)])
  df$isEnd <- c(!df$isConsec[-nrow(df)], T)
  bed <- data.frame(start=df$pos[df$isStart], end=df$pos[df$isEnd])
  bed$chr <- myChr
  return(bed[, c("chr", "start", "end")])
}

# From the umi4cPackage
# https://bitbucket.org/tanaylab/umi4cpackage/src/59c73399f73d0b6d14dd03c18317ab4d0569b332/R/p4c_package.R#lines-1015
# Calc geometric mean of coordinates relative to distance from bait
.p4cWinGeoMeanCoordinate <- function(coords, scale, bait_x){
  bait_idx <- which(abs(coords - as.numeric(bait_x)) == min(abs(coords - as.numeric(bait_x))))[1]
  offsets <- abs(coords - coords[bait_idx])
  offsets[bait_idx] <- 1  #avoiding plugging 0s in the log
  
  
  mean_offsets1 <- exp(zoo::rollmean(log(offsets[1:bait_idx]), scale, fill = NA))
  mean_offsets2 <- exp(zoo::rollmean(log(offsets[(bait_idx + 1):length(offsets)]), scale, 
                                     fill = NA))
  
  # deal with NAs near the bait
  near_bait1 <- sapply((bait_idx - scale/2 + 1):bait_idx, function(i) exp(mean(log(offsets[i:bait_idx]))))
  near_bait2 <- sapply((bait_idx + 1):(bait_idx + ceiling(scale/2) - 1), 
                       function(i) exp(mean(log(offsets[(bait_idx + 1):i])))) #If even scale it is not working
  
  mean_offsets1[(bait_idx - scale/2 + 1):bait_idx] <- near_bait1
  mean_offsets2[1:(ceiling(scale/2) - 1)] <- near_bait2 #If even scale it is not working
  
  mean_offsets1 <- mean_offsets1 * -1
  mean_offsets1[bait_idx] <- 0
  
  mean_offsets <- c(mean_offsets1, mean_offsets2)
  m_coords <- coords[bait_idx] + mean_offsets
  
  return(m_coords)
}

# Inspired from umi4cpackage
# https://bitbucket.org/tanaylab/umi4cpackage/src/59c73399f73d0b6d14dd03c18317ab4d0569b332/R/p4c_package.R#lines-779:783
rollMean <- function(prof, scale){
  na.pos <- which(is.na(prof))
  prof[na.pos] <- 0
  has.value <- rep(1, length(prof))
  has.value[na.pos] <- 0
  cprof <- cumsum(prof)
  cprof.values <- cumsum(has.value)
  if (scale %% 2 == 0){
    d <- scale / 2
    sums <- c(rep(NA, d), 
              (cprof[(2 * d + 1):length(cprof)] - cprof[1:(length(cprof) - 2 * d)]),
              rep(NA, d))
    nb.values <- c(rep(NA, d), 
                   (cprof.values[(2 * d + 1):length(cprof.values)] - cprof.values[1:(length(cprof.values) - 2 * d)]),
                   rep(NA, d))
  } else {
    d <- (scale - 1) / 2
    sums <- c(rep(NA, d),cprof[scale]/scale, 
              (cprof[(scale + 1):length(cprof)] - cprof[1:(length(cprof) - scale)])/scale,
              rep(NA, d))
    nb.values <- c(rep(NA, d), cprof.values[scale] / cprof.values[scale], 
                   (cprof.values[(scale + 1):length(cprof.values)] - cprof.values[1:(length(cprof.values) - scale)])/scale,
                   rep(NA, d))
  }
  return(sums / nb.values)
}
