# Plot idea from http://dx.doi.org/10.1016/j.cell.2014.06.049
# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)

output <- read.table("tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
withEdges <- FALSE # Bernoulli vs. Graphical
alpha <- 0.5 # Alpha 0 or 0.5

# Load cluster assignment
if (withEdges) {
    load(file = "clusteringEMk22alpha0.5best.RData")
    newclustermembershipbest <- newclustermembership
} else if (alpha == 0.5){
    load('clusteringEMemptymutk22alpha0.5tau340.RData')
} else if (alpha == 0) {
    load('clusteringEMemptymutk22alpha0.001tau340.RData')
} 

load('../../TCGAnewjk.rData') # Load cancer types

# Map clinical data to clustering
groupId <- match(output$id, combyordered$id)
clinical <- cbind(output, newclustermembershipbest[groupId], combyordered$type[groupId])
colnames(clinical) <- c(colnames(output), 'group', 'tissue')
# Change group from 1:22 to A:V
clinical$group <- as.factor(clinical$group)
levels(clinical$group) <- LETTERS[1:length(unique(newclustermembershipbest))]
clinical <- subset(clinical, !is.na(clinical$group))
table(clinical$group)
clinical[,c('age','event','time')] <- lapply(clinical[,c('age','event','time')], as.numeric)
# Remove subtype from tissue
#levels(clinical$tissue) <- sapply(levels(clinical$tissue), function(tissue) {strsplit(tissue, " ")[[1]][1]})
clinical[,c('stage','t','n','m','group','tissue')] <- lapply(clinical[,c('stage','t','n','m','group','tissue')], factor)

##########################################
# Change in likelyhood ratio test when variable are added (stage)
##########################################

# Clinical only
stageTest <- summary(coxph(Surv(time, event) ~ stage, data = clinical, na.action = "na.omit"))$logtest[1]
ageTest <- summary(coxph(Surv(time, event) ~ stage + age, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue
tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue + group
groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + group, data = clinical, na.action = "na.omit"))$logtest[1]


mat <- matrix(c(stageTest, ageTest-stageTest, tissueTest-ageTest, groupTest-tissueTest), nrow = 4)
mat <- mat/2 # Log-rank statistic is twice the chisquare statistic
col <- rev(brewer.pal(4, "Greys"))
if (withEdges) {
    pdf('cluster_improved_survival_prediction_k_22_alpha_0.5.pdf', width=5, height = 10)
} else if (alpha == 0){
    pdf('cluster_improved_survival_prediction_k_22_alpha_0.001_empty', width=5, height = 10)
} else if (alpha == 0.5) {
    pdf('cluster_improved_survival_prediction_k_22_alpha_0.5_empty', width=5, height = 10)
} 
barplot(mat,col = col, ylim=c(0,ceiling(sum(mat)/100)*100), ylab='Change in LR statistic')
# P-value chisquare distribution
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE), 10)
legend("center", legend = rev(c('Clinical (stage)','Clinical (age)','Tissue type',paste0('Cluster (LR=',LR,', p=',pvalue,')'))),
       fill = rev(col), bg = 'white', cex = 0.7)
dev.off()


##########################
# Test specific groups

coxph.subset <- function(type, g1, g2, data) {
    clinical.subset <- subset(x = data, subset = tissue == type & group %in% c(g1,g2))
    clinical.subset$group <- factor(clinical.subset$group)
    summary(coxph(Surv(time, event) ~ group, data = clinical.subset, na.action = "na.omit"))
}

# G vs. N in lgg
coxph.subset("lgg","N","G",clinical)

# I vs. K in caodread
coxph.subset("coadread","I","K",clinical)

# L vs. V in laml
coxph.subset("laml", "L", "V", clinical)

# U vs. V in ov
coxph.subset("ov", "U", "V", clinical)



#########################
# Without correction for confounding factors
#########################
clusterTest <- summary(coxph(Surv(time, event) ~ group, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- clusterTest/2
pvalue <- pchisq(q = clusterTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE)

paste(round(LR, 1), "&", pvalue, sep = " ")
