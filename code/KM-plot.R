# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)
output <- read.table("../data/tcga-clinical-information.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Load classification by mutation profile (cluster assignment)
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

# Kaplan-Meier curve for 22 groups
# Colors A -- V
colourysdots<-c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
                    "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
                    "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
                    "#77AADD","#771155","#AA4488","#CC99BB")

pdf('./overall-survival-cluster-types-22.pdf', width=10, height=10)
par(mar = c(4.5,4.5,0.5,0.5))
    plot(survfit(Surv(time = as.numeric(time)/365, event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', cex.lab=1.75, cex.axis=1.5, mark.time = T, xlim=c(0,24))
  legend(x = 19, y = 1.01, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1.75, ncol = 2,title="Cluster", bg='white')
    
    xs<-c(15,5.6,10,15,9.2,11.5,15.2,18.5,11.5,10.9,8.9,9.3,9.6,13.7,14.4,16.5,15,12.6,11.1,20.8,22.5,21.5)
    ys<-c(0.86,0.27,0.34,0.645,0.385,0.445,0.115,0.125,0.505,0.03,0.425,0.2,0.545,0.1,0.28,0.17,0.4,0.03,0.32,0.28,0.125,0.41)
    for(ii in 1:length(xs)){
      text(xs[ii],ys[ii],col=colourysdots[ii],LETTERS[ii],cex=1.5)
    }
    
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
    print(xtable(cbind(testStat$coefficients[,2], testStat$conf.int[,c(3:4)], testStat$coefficients[,5]), digits = c(0,3,3,3,-1)))
}


xtPrint(clinical, "V") # Baseline V (biggest group as baseline)

# Contributions by cancer type for cluster E
#boxplot(as.numeric(clinical$time[clinical$group == 'E']) ~ as.factor(clinical$type[clinical$group == 'E']), las=2)

###################################################################
# Cox proportional hazards with different baselines, stratified for stage, age and type
###################################################################
xtPrint <- function(clinical, ref) {
    clinical$group <- relevel(clinical$group, ref=ref)
    testStat <- summary(coxph(Surv(as.numeric(time)/365, as.numeric(event)) ~ stage + as.numeric(age) + tissue + group, data = clinical)) # Baseline A
    testStat$coefficients[,"Pr(>|z|)"] <- p.adjust(testStat$coefficients[,"Pr(>|z|)"], method = "fdr")
    print(xtable(cbind(testStat$coefficients[,2], testStat$conf.int[,c(3:4)], testStat$coefficients[,5]), digits = c(0,3,3,3,-1)))
}
xtPrint(clinical, "V") # Baseline V


###################################################################
# Comparing specific groups (laml E & J, brca J & K, or ucec a & D)
###################################################################
coxph.subset <- function(tissuetoselect, g1, g2, data) {
    clinical.subset <- subset(x = data, tissue == tissuetoselect & group %in% c(g1,g2))
    clinical.subset$group <- factor(clinical.subset$group)
    summary(coxph(Surv(time, as.numeric(event)) ~ group, data = clinical.subset, na.action = "na.omit"))
}
# G vs. N in lgg
coxph.subset("lgg","N","G",clinical)

# I vs. K in coadread
coxph.subset("coadread","I","K",clinical)

# L vs. V in laml
coxph.subset("laml", "L", "V", clinical)

