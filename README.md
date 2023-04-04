# Tumor_Heterogeneity

This repository contains the R scripts calculating the tumor heterogeneity and applying patients clustering introduced in paper “Evolutionary states and trajectories characterized by distinct pathways stratify ovarian high-grade serous carcinoma patients”.

1. The "cellular_freqs.R" calculates the intra and inter-tumor heterogeneity using clonal mutational frequencies as input.

2. The "kmeansw.R" script does the clustering based on pre-assigned centroids described in the Method section of the paper.

3. Regarding the clustering in validation set, meanig assigning new patients to old clusters, we do the clustering with the discovery data *only*. This is done in "test2-new.R". One will need to grab the clustering variable "clu". This will have the labels clu$L (for the discovery data) and the
cluster parameters in clu$Y.

4. We the beginning of "test2-new.R" to compute "X" and "W" (that is, avg$hcp and avg$dup) for the validation data (can include the discovery data too).

5. Using the old "clu$Y" and the new "X", "W" we map the two together. This is done in "kmeansw.R" L61-70. Now, new.clu$L will contain the cluster labels for the new data. 


## Contributors
Sanaz Jamalzadeh   
Antti Häkkinen
