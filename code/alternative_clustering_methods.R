# Compare clustering method with alternative methods

rm(list = ls())
mut.matrix <- read.table("../data/binary-mutation-matrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


###########################
# Hierarchical clustering
###########################
hamming <- function(X) {
    D <- (1 - X) %*% t(X)
    D + t(D)
}

h.dist <- hamming(data.matrix(mut.matrix))
dist <- as.dist(h.dist)
clusters <- hclust(dist)
clustermembership.hclust <- cutree(clusters, k = 22)

###########################
# NMF clustering (remove zero lines and assign to specific group for survival analysis)
###########################
library(NMF)
mat <- t(mut.matrix)
mat <- mat[,colSums(mat) > 0]
set.seed(1234)
result <- nmf(x = mat, rank=21)
nmf_clusters <- apply(coef(result), 2, which.max)
clustermembership.nmf <- nmf_clusters

###########################
# K-Means
###########################
set.seed(1234)
# Best of 240 iterations
result <- kmeans(x = mut.matrix, centers = 22, nstart = 240)
clustermembership.kmeans <- fitted(result, method = "classes")

###########################
# mclust
###########################
library(mclust)
BIC <- mclustBIC(mut.matrix, G = 22)
result <- Mclust(mut.matrix, 22, x = BIC)
clustermembership.mclust <- summary(result, parameters = TRUE)$classification

###########################
# clean-up and save
###########################
rm(hamming, h.dist, dist, clusters, mat, result, nmf_clusters, BIC)
save.image("other_clustering.Rdata")
