
# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)

output <- read.table("../data/tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cluster.membership <- read.table("../data/annotation-matrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cluster.membership <- cluster.membership[,c(1:4)] # select columns

#rename some columns
colnames(cluster.membership)[colnames(cluster.membership)=="type"] <- "tissue"
colnames(cluster.membership)[colnames(cluster.membership)=="cluster"] <- "group"

# Change group from 1:22 to A:V
cluster.membership$group <- as.factor(cluster.membership$group)
levels(cluster.membership$group) <- LETTERS[1:length(unique(cluster.membership$group))]

# Merge clinical data with clustering
clinical <- merge(output,cluster.membership, by="id", sort=FALSE)

table(clinical$group)

clinical[,c('age','event','time')] <- lapply(clinical[,c('age','event','time')], as.numeric)
# Remove subtype from tissue
clinical[,c('stage','group','tissue')] <- lapply(clinical[,c('stage','group','tissue')], factor)

##########################################
# Change in likelihood ratio test when variables are added
##########################################

# Clinical only
stageTest <- summary(coxph(Surv(time, event) ~ stage, data = clinical, na.action = "na.omit"))$logtest[1]
ageTest <- summary(coxph(Surv(time, event) ~ stage + age, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue
tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue, data = clinical, na.action = "na.omit"))$logtest[1]
# Clinical + tissue + group
groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ stage + age + tissue + group, data = clinical, na.action = "na.omit"))$logtest[1]


LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE), 10)

paste(round(LR, 1), "&", pvalue, sep = " ")


#########################
# Without correction for confounding factors
#########################

clusterTest <- summary(coxph(Surv(time, event) ~ group, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- clusterTest/2
pvalue <- pchisq(q = clusterTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE)

paste(round(LR, 1), "&", pvalue, sep = " ")
