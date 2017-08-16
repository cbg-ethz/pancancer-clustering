# Compare clustering method with alternative methods

rm(list = ls())
load("TCGAnewjk.rData")


###########################
# Hierarchical clustering
###########################
hamming <- function(X) {
    D <- (1 - X) %*% t(X)
    D + t(D)
}

h.dist <- hamming(data.matrix(combyordered[,-1:-2]))
dist <- as.dist(h.dist)
clusters <- hclust(dist)
clustermembership.hclust <- cutree(clusters, k = 22)
table(clustermembership.hclust)

###########################
# NMF clustering (remove zero lines and assign to specific group for survival analysis)
###########################
library(NMF)
mat <- t(combyordered[,-1:-2])
mat <- mat[,colSums(mat) > 0]
set.seed(1234)
result <- nmf(x = mat, rank=21)
nmf_clusters <- apply(coef(result), 2, which.max)
clustermembership.nmf <- nmf_clusters
table(clustermembership.nmf)

###########################
# K-Means
###########################
set.seed(1234)
# Best of 240 iterations
result <- kmeans(x = combyordered[,-1:-2], centers = 22, nstart = 240)
clustermembership.kmeans <- fitted(result, method = "classes")

###########################
# mclust
###########################
library(mclust)
BIC <- mclustBIC(combyordered[,-1:-2], G = 22)
result <- Mclust(combyordered[,-1:-2], 22, x = BIC)
clustermembership.mclust <- summary(result, parameters = TRUE)$classification

###########################
# clean-up and save
###########################
rm(combyordered, hamming, h.dist, dist, clusters, mat, result, nmf_clusters, BIC)
save.image("other_clustering.Rdata")
