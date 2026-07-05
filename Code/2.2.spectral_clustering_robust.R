### Running spectral clustering on the CT data
## This script is run after Running_parameter_iterations.r in order to determine the optimal 
# cluster number, hyperparameter (alpha), and nearest neighbors parameter (K)

library(SNFtool)
library(cluster)
library(MASS)
library(ggplot2)
library(corrplot)
library(dplyr)
library(broom)
library(tidyr)
library(psych)
library(dunn.test)
library(fossil)
setwd("~/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/clustering")
# setting output directory
directory = ("~/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/clustering/ABCD")

# source code for function to determine data integration and clustering across resampling using SNF and spectral clustering
source("./function_robust_core-clustering.R")

# set up for clustering analysis
# importing the different data types as individual participant matrices 
data = read.csv("./ABCD/ct_z_adhd.csv",header = TRUE)
data_adhd = filter(data, group == "ADHD")
subjects  = data.frame(id = data_adhd$id) #list of participants
brain = data_adhd[,c(8:226)] # cortical thickness (n=219)
# normalizing measures within each data type using a function from the SNF package
brain = standardNormalization(brain)

# setting the parameters (finalized after comparisons using Running_parameter_iterations.r )
K =30;		# number of neighbors, usually (10~30), usually sample size/10
alpha = 0.8;  	# hyperparameter, usually (0.3~0.8)

# creating participant distance matrices using euclidean distances
Dist_brain = dist2(as.matrix(brain),as.matrix(brain))
write.csv(Dist_brain,file=file.path(directory, paste("ct_dist_ABCD.csv", sep="")))

# creating participant affinity matrices within each data type
AM_brain = affinityMatrix(Dist_brain,K,alpha) 

### Calculating similarity matrix and spectral clustering groups
# Clustering
C = 2 #setting cluster number
subgroup = spectralClustering(AM_brain,C)
displayClustersWithHeatmap(AM_brain, subgroup)
# calling function to cluster participants using spectral clustering across resampling 80% of participants 1000 times
robust.W = RobustCoreClusteringMatrix(feature.affinity.mat.list = list(AM_brain),exp.num.samples = 1000, num.clusts = C)
#Two matrices - Dense Core Cluster Matrix and Sparse Core Cluster Matrix
dense <- robust.W[1]
dense <- matrix(unlist(dense), ncol = 1112, byrow = TRUE)
sparse <- robust.W[2]
sparse <- matrix(unlist(sparse), ncol = 1112, byrow = TRUE)

# displaying clusters
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))

# saving an image of the cluster heatmap
png('./dense_0.8_30.png',width = 1500, height = 1500, units = "px", bg = "white", res = 300)
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
dev.off()

png('./sparse_0.8_30.png',width = 1500, height = 1500, units = "px", bg = "white", res = 300)
displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))
dev.off()

#saving clustering similarity matrix
write.matrix(dense, file=file.path(directory, paste("dense_k30_0.8_matrix.csv", sep="")),sep = ",")
write.matrix(sparse, file=file.path(directory, paste("sparse_k30_0.8_matrix.csv", sep="")),sep = ",")

# Find cluster labels of individuals using the robust clustering similarity matrix
robust.groups.df = RobustCoreClusteringClusters(core.clustering.list = robust.W,num.clusts = C,verbose = TRUE)
clusters <- cbind(subjects, robust.groups.df$groups)
colnames(clusters)[2] = 'subtype'
table(clusters$subtype)

## calculating silouette width for each participant and the silhouette plot
dissim <- 1 - dense
dissim <- as.matrix(dissim)
clusters$subtype <- as.integer(clusters$subtype)
sil <- silhouette(clusters$subtype, dmatrix = dissim)
summary(sil)

# saving the silhouette plot
plot(sil, col = c("palevioletred", "sandybrown"))
tiff("Silhouette_CT.tiff", width = 18, height = 18, units = "cm", bg="white",res = 200)
plot(sil, col = c("palevioletred", "sandybrown"))
dev.off()

# adding silhouette widths for each individual
clusters$silhouette_width <- 0
for (i in 1:1112){
  clusters[i, 3] <- sil[i ,3]
}
write.csv(clusters, file=file.path(directory, paste("2clust_groups_CT_ABCD.csv", sep="")))




