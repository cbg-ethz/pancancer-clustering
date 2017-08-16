rm(list = ls())
library(limma)
library(biomaRt)
mut.matrix <- read.table("../data/binary-mutation-matrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
genes <- colnames(mut.matrix)
# Downloaded from https://string-db.org/ (Homo sapiens only)
network <- read.table("9606.protein.links.v10.5.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
network[,1:2] <- data.frame(cbind(substring(network$protein1, 6), substring(network$protein2, 6)), stringsAsFactors=FALSE)

####################
# Map Ensemble protein to symbols
ensbl.network <- unique(c(network[,1],network[,2]))
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id","hgnc_symbol"),
                 filters = "ensembl_peptide_id", values = ensbl.network, mart = mart)

####################
# Only keep edges with genes in mutation matrix
network.genes <- data.frame(cbind(results$hgnc_symbol[match(network$protein1, results$ensembl_peptide_id)], 
                                  results$hgnc_symbol[match(network$protein2, results$ensembl_peptide_id)]), stringsAsFactors = FALSE)
network.genes.clean <- network.genes[!is.na(network.genes[,1]) & !is.na(network.genes[,2]),]
symbols <- alias2SymbolTable(genes)[!is.na(alias2SymbolTable(genes))]
network.clean <- network.genes.clean[alias2SymbolTable(network.genes.clean[,1]) %in% symbols & alias2SymbolTable(network.genes.clean[,2]) %in% symbols,]

rm(network, ensbl.network, results, network.genes, network.genes.clean, mart)

####################
# Adjacency matrix of string network
network.clean.adj <- data.frame(matrix(data = 0, nrow = length(symbols), ncol = length(symbols), dimnames = list(symbols, symbols)), stringsAsFactors = FALSE)
invisible(apply(network.clean, 1, function(tuple) { network.clean.adj[tuple[1], gsub("-","\\.",tuple[2])] <<- 1 })) # colname changes - to .

####################
# Create adjacency matrix from 
id.missing.gene <- which(is.na(alias2SymbolTable(genes)))
types <- c("blca","brca","cesc","coadread","esca","gbm","hnsc","kirc","kirp","laml","lgg","lihc","luad","lusc","ov","paad","pcpg","prad","sarc","stad","thca","ucec")
cancer.adj <- lapply(types, function(type,id) { 
    mat <- read.table(paste("../data/networks/", type, "-network.txt", sep = ""), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    return(mat[-id,-id]) 
    }, id.missing.gene)

###################
# Use sum of all cancers
cancer.adj.all <- 1*(Reduce('+', cancer.adj) > 0) # Sum up adj matrices of all cancer types
cancer.adj.all <- cancer.adj.all | t(cancer.adj.all) # Make undirected

permutation.test <- function(BNmatrix,interactionmatrix) {
    permy <- sample(200)
    permymatrix<-interactionmatrix[permy,permy]
    mn <- mean(BNmatrix==permymatrix)
    return(mn)
}

set.seed(1234)
n.rep <- 1e06
mean.sim.network <- replicate(n.rep, permutation.test(cancer.adj.all, network.clean.adj))
result <- sum(mean.sim.network > mean(cancer.adj.all == network.clean.adj))/n.rep
save.image("result.rdata")
