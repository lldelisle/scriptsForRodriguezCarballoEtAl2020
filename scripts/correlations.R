## Load needed packages
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("gridExtra")
safelyLoadAPackageInCRANorBioconductor("GGally")

######## Functions #########

getFullMatrix <- function(values){
  # Get the number of bins
  n <- (1 + sqrt(length(values) * 8 + 1)) / 2
  # Make the skeleton:
  df <- expand.grid(bin1 = 0:(n - 1), bin2 = 0:(n - 1), stringsAsFactors = F)
  df$value <- 0
  df <- df[order(df$bin1, df$bin2), ]
  df[df$bin1 < df$bin2, "value"] <- values
  combin <- paste(df$bin1, df$bin2, sep = "_")
  rownames(df) <- combin
  combinT <- paste(df$bin2, df$bin1, sep = "_")
  df$value <- df$value + df[combinT, "value"]
  return(df)
}

########

plotDirectory <- commandArgs(TRUE)[1]

# Get all file names
all.vhic <- system("ls */model_*_onTDOM/*/*/HoxD_final_output_*_*_*/vhic_HoxD.txt", intern=TRUE)

# I assume the files are:
# *_[exclusion_size][iteration]/model_[experiment_name]_onTDOM/*/*/HoxD_final_output_[uZ]_[lZ]_[max.d]/vhic_HoxD.txt

# First extract the condition of input
exclusion.iteration <- sapply(all.vhic, function(s){strsplit(s, "/")[[1]][1]})
iteration <- sapply(exclusion.iteration, function(s){substring(s, nchar(s))})
exclusion <- sapply(exclusion.iteration, function(s){as.numeric(substring(strsplit(s, "_")[[1]][2], 1, 1))})
experiment <- sapply(all.vhic, function(s){gsub("model_", "", gsub("_onTDOM", "", strsplit(s, "/")[[1]][2]))})

# Then get the output properties
finalFolder <- sapply(all.vhic, function(s){gsub("HoxD_final_output_", "", strsplit(s, "/")[[1]][5])})
properties <- matrix(as.numeric(unlist(strsplit(finalFolder, "_"))), ncol = 3, byrow = T)
colnames(properties) <- c("upper.Z", "lower.Z", "max.distance")

# Combine everything
all.properties <- data.frame(file=all.vhic, i=iteration, exclusion=exclusion, experiment=experiment, properties)

# Adjust the names
all.properties$expe.exclu <- paste(all.properties$experiment, all.properties$exclusion, sep = "_")
all.properties$expe.exclu.i <- paste(all.properties$experiment, all.properties$exclusion, all.properties$i, sep = "_")
rownames(all.properties) <- all.properties$expe.exclu.i
names(all.vhic) <- all.properties$expe.exclu.i

# Reorder
all.vhic <- all.vhic[sort(names(all.vhic))]
all.properties <- all.properties[sort(names(all.vhic)), ]

# Display some outputs
pdf(paste0(plotDirectory, "/output_boundaries.pdf"), title = "output boundaries")
g <- ggplot(all.properties, aes(upper.Z, i)) +
  geom_point() + 
  facet_grid(rows = vars(experiment), cols = vars(exclusion)) +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle(label = "Upper Z")
print(g)
g <- ggplot(all.properties, aes(lower.Z, i)) +
  geom_point() + 
  facet_grid(rows = vars(experiment), cols = vars(exclusion)) +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle(label = "Lower Z")
print(g)
g <- ggplot(all.properties, aes(max.distance, i)) +
  geom_point() + 
  facet_grid(rows = vars(experiment), cols = vars(exclusion)) +
  theme(strip.text.y = element_text(angle = 0)) +
  ggtitle(label = "Maximum distance")
print(g)
dev.off()

# Read all files
all.vhic.m <- lapply(all.vhic, read.delim, h=F, sep=",", col.names=c("bin1", "bin2", "val"))

# Simplify
all.vhic.m <- lapply(all.vhic.m, subset, subset=bin1 < bin2)
all.vhic.v <- lapply(all.vhic.m, function(df){df[, 3]})

# Group per experiment
all.vhic.per.expe <- lapply(unique(all.properties$experiment),
                            function(e){
                              do.call(cbind, subsetByNamesOrIndices(all.vhic.v, rownames(subset(all.properties, experiment == e))))
                            })
names(all.vhic.per.expe) <- unique(all.properties$experiment)

# Check how many iteration I have for each experiment
sapply(all.vhic.per.expe, ncol)

# Perform correlation within each experiment
correlations <- lapply(all.vhic.per.expe, cor, method="spearman")

# In order to display the same scale, get the minimum
min.value <- min(unlist(correlations))

# The clusters and the mean correlation between them will be stored
all.properties$cluster <- factor(1, levels = 1:nrow(all.properties))
cluster.distances <- NULL

# Will make a multi plot
plot_list=list()
for (i in 1:length(correlations)){
  # First perform a hierarchical clustering
  hc <- hclust(dist(correlations[[i]]), method = "ward.D2")
  k <- 1
  max.means <- 0
  # Cut the clustering to get 1 cluster
  clusters <- cutree(hc, k)
  all.pairs <- NULL
  # I will increase the number of clusters until
  # the mean correlation between 2 clusters is above 0.95
  while(max.means < 0.95){
    # Store values of previous step
    previous.clusters <- clusters
    previous.all.pairs <- all.pairs
    # Increase the number of cluster
    k <- k + 1
    # Determine the clusters
    clusters <- cutree(hc, k)
    # Get all possible pairs of cluster
    all.pairs <- as.data.frame(t(combn(k, 2)))
    # evaluate the mean of correlation
    all.pairs$means <- apply(all.pairs, 1, function(v){mean(correlations[[i]][names(clusters[clusters == v[1]]),
                                                                              names(clusters[clusters == v[2]])])})
    # Get the maximum value
    max.means <- max(all.pairs$means)
  }
  # We need to go back one step
  n.clusters <- k - 1
  clusters <- previous.clusters
  # We add the info to all.properties and cluster.distances
  all.properties[names(clusters), "cluster"] <- clusters
  if (!is.null(previous.all.pairs)){
    previous.all.pairs$experiment <- names(correlations)[i]
    cluster.distances <- rbind(cluster.distances, previous.all.pairs)
  }
  # We plot a pheatmap showing both the clusters and the exclusion size
  p <- pheatmap(correlations[[i]],
                annotation_row = subset(all.properties, select=exclusion),
                annotation_col = all.properties[, c("exclusion", "cluster")],
                show_rownames = F,
                show_colnames = F,
                main = names(correlations)[i],
                breaks = seq(min.value, 1, length.out = 100),
                annotation_names_row = F,
                annotation_names_col = F,
                clustering_method = "ward.D2",
                annotation_legend = T,
                silent = T
  )
  # We add the plot to the plot_list
  plot_list[[i]] <- p[[4]]
}
# We arrange the plot_list
g <- grid.arrange(arrangeGrob(grobs=plot_list, ncol=sqrt(length(correlations))))
# Save it in the plotDirectory
ggsave(paste0(plotDirectory, "/correlations_wd2.pdf"), g, width = 17, height = 17)

# Now we will perform averages:

# Do mean total:
mean.v <- lapply(all.vhic.per.expe, rowMeans)

# Do mean per cluster and write all averages:
for (e in unique(all.properties$experiment)){
  ## Total:
  # We need to rebuild the full matrix
  df <- getFullMatrix(mean.v[[e]])
  # And write it as csv without colnames
  write.table(x = df, file = paste0("average__", e, "__all.txt"), row.names = F, col.names = F, sep = ",")
  ## Per cluster:
  nb.clusters <- max(as.numeric(as.character(all.properties[all.properties$experiment == e, "cluster"])))
  # We only write it if there is more than 1 cluster
  if (nb.clusters > 1){
    for (k in 1:nb.clusters){
      iterations.in.cluster <- rownames(subset(all.properties, experiment == e & cluster == k))
      # We perform the mean on values
      if (length(iterations.in.cluster) > 1){
        values <- rowMeans(all.vhic.per.expe[[e]][, iterations.in.cluster])
      } else {
        values <- all.vhic.per.expe[[e]][, iterations.in.cluster]
      }
      # Rebuild the full matrix
      df <- getFullMatrix(values)
      # Write it
      write.table(x = df, file = paste0("average__", e, "__cluster", k, ".txt"), row.names = F, col.names = F, sep = ",")
    }
  }
}

# Output the all.properties
write.table(all.properties, paste0(plotDirectory, "/4Cin_all_properties.txt"),
            sep = "\t", quote = F, row.names = F)

# Output the cluster distances
colnames(cluster.distances) <- c("cluster1", "cluster2", "mean", "experiment")
write.table(cluster.distances, paste0(plotDirectory, "/4Cin_cluster_distances.txt"),
            sep = "\t", quote = F, row.names = F)
