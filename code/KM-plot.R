# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)
library(gplots)
output <- read.table("../data/tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load classification by mutation profile (cluster assignment)
cluster.membership <- read.table("../data/annotation-matrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Map clinical data to clustering
groupId <- match(output$id, rownames(cluster.membership))
clinical <- cbind(output, cluster.membership$cluster[groupId], cluster.membership$type[groupId])
colnames(clinical) <- c(colnames(output),'group','type')
table(clinical$group)
clinical$group <- factor(as.numeric(clinical$group))
levels(clinical$group) <- LETTERS[seq(length(levels(clinical$group)))]
# Remove NAs
clinical <- subset(clinical, !is.na(clinical$group) & !is.na(clinical$event))
clinical$group <- factor(clinical$group)
table(clinical$group)
levels(clinical$type) <- sapply(levels(clinical$type), function(tissue) {strsplit(tissue, " ")[[1]][1]})

# Kaplan-Meier curve for 22 groups
# Colors A -- V
colourysdots<-c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
                    "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
                    "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
                    "#77AADD","#771155","#AA4488","#CC99BB")
pdf('../overall-survival-cluster-types-22.pdf', width=10, height=10)
par(mar = c(4.5,4.5,0.5,0.5))
    plot(survfit(Surv(time = as.numeric(time)/365, event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', cex.lab=1.5, cex.axis=1.5, mark.time = T)
    legend(x = 21.5, y = 1.01, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1.75, ncol = 2,title="Cluster", bg='gray90')
dev.off()

###################################################################
# 5-year overall survival
###################################################################
clinical.V <- subset(x = clinical, subset = clinical$group == "V")
fit <- survfit(Surv(time = as.numeric(time)/365, event = as.numeric(event)) ~ group, data = clinical)
summary(fit, times = 5)


###################################################################
# Cox proportional hazards with different baselines
###################################################################
library(xtable)
xtPrint <- function(clinical, ref) {
    clinical$group <- relevel(clinical$group, ref=ref)
    testStat <- summary(coxph(Surv(as.numeric(time)/365, as.numeric(event)) ~ group, data = clinical))
    testStat$coefficients[,"Pr(>|z|)"] <- p.adjust(testStat$coefficients[,"Pr(>|z|)"], method = "fdr")
    print(xtable(cbind(testStat$coefficients[,2], testStat$conf.int[,c(3:4)], testStat$coefficients[,5]), digits = c(0,3,3,3,-2)))
}


xtPrint(clinical, "V") # Baseline V (biggest group as baseline)

# Contributions by cancer type for cluster E
#boxplot(as.numeric(clinical$time[clinical$group == 'E']) ~ as.factor(clinical$type[clinical$group == 'E']), las=2)

###################################################################
# Cox proportional hazards with different baselines, stratified for stage, age and type
###################################################################
xtPrint <- function(clinical, ref) {
    clinical$group <- relevel(clinical$group, ref=ref)
    testStat <- summary(coxph(Surv(as.numeric(time)/365, as.numeric(event)) ~ stage + as.numeric(age) + type + group, data = clinical)) # Baseline A
    testStat$coefficients[,"Pr(>|z|)"] <- p.adjust(testStat$coefficients[,"Pr(>|z|)"], method = "fdr")
    print(xtable(cbind(rep(NA, length(testStat$coefficients[,2])), testStat$coefficients[,2], testStat$conf.int[,c(3:4)], testStat$coefficients[,5]), digits = c(0,0,3,3,3,-2)))
}
xtPrint(clinical, "V") # Baseline V


###################################################################
# Comparing specific groups (laml E & J, brca J & K, or ucec a & D)
###################################################################
coxph.subset <- function(tissue, g1, g2, data) {
    clinical.subset <- subset(x = data, subset = tissue == type & group %in% c(g1,g2))
    clinical.subset$group <- factor(clinical.subset$group)
    summary(coxph(Surv(time, as.numeric(event)) ~ group, data = clinical.subset, na.action = "na.omit"))
}
# G vs. N in lgg
coxph.subset("lgg","N","G",clinical)

# I vs. K in caodread
coxph.subset("coadread","I","K",clinical)

# L vs. V in laml
coxph.subset("laml", "L", "V", clinical)

# U vs. V in ov
coxph.subset("ov", "U", "V", clinical)
