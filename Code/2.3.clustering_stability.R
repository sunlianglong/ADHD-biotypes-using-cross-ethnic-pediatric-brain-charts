### Script to run resampling stability analysis

library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

setwd('~/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/clustering')
source("./function_bootstrapping.R")

data = read.csv("./ct_z_merged.csv",header = TRUE)
data_adhd = filter(data, group == "ADHD")
subjects  = data.frame(id = data_adhd$id) #list of participants
CT = data_adhd[,c(6:224)]
directory <- ("~/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/clustering/clustering_stability/")

numboot=1000
nsub=587
K=30
alpha=0.8
bootsize=0.8
clusters=2

# Setting up which participants will be included in each permutation
permutation_matrix <- bootstrapping_SNF(numboot=numboot, nsub=nsub, bootsize=bootsize)
# Getting clustering solutions for all the permutations of sampled participants using SNF 
clus_sil <- clustering(perms=permutation_matrix, bootsize=bootsize, K=K, alpha=alpha, clusters=clusters, CT=CT)
# Dividing output matrix into clusters and silhouette widths
clus_out <- clus_sil[1:(numboot), ]
silhouette_width <- clus_sil[(numboot+1):(numboot*2), ]
# getting the adjusted rand index between all clustering solutions
list_randindex <- stability(clus_out=clus_out, perms=permutation_matrix) # returns b rand index:adjusted rand index
list_adjustedrandindex <- list_randindex[ ,1001:2000]
# Calculate how often each participant is clustered together and the probability that they will be clustered together
percent_agree <- percent_agree(clus_out=clus_out)

write.csv(permutation_matrix, file=file.path(directory, paste("permutation_matrix_1000perms.csv", sep="")))
write.csv(percent_agree, file=file.path(directory, paste("Percent_agree_1000perms.csv", sep="")))
write.csv(list_adjustedrandindex, file=file.path(directory, paste("adjrandindex_1000perms.csv", sep="")))
write.csv(list_randindex, file=file.path(directory, paste("Rand_indices_1000perms.csv", sep="")))
write.csv(clus_out, file=file.path(directory, paste("clu_solution_1000perms.csv", sep="")))
write.csv(silhouette_width, file=file.path(directory, paste("silhouette_1000perms.csv", sep="")))

