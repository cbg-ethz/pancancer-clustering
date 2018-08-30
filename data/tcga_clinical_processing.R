
clinical.data.all <- read.csv("tcga_clinical_Cell2018.csv")

#rename some columns
colnames(clinical.data.all)[colnames(clinical.data.all)=="bcr_patient_barcode"] <- "id"
colnames(clinical.data.all)[colnames(clinical.data.all)=="age_at_initial_pathologic_diagnosis"] <- "age"
colnames(clinical.data.all)[colnames(clinical.data.all)=="ajcc_pathologic_tumor_stage"] <- "stage"
colnames(clinical.data.all)[colnames(clinical.data.all)=="OS"] <- "event"
colnames(clinical.data.all)[colnames(clinical.data.all)=="OS.time"] <- "time"

clinical.data.selected <- clinical.data.all[,c("id", "age", "stage", "event", "time")]

clinical.data.selected <- clinical.data.selected[-which(clinical.data.selected$time=="#N/A"),]
clinical.data.selected <- clinical.data.selected[-which(clinical.data.selected$age=="#N/A"),]

# select only those patients in our cohort

cluster.membership <- read.table("annotation-matrix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cluster.membership.id <- data.frame(id=cluster.membership$id)

clinical.data <- merge(clinical.data.selected, cluster.membership.id, by="id", sort=FALSE)

library(plyr)

names1 <- paste0("Stage I",c("","A", "B"))
names2 <- rep("stage i", length(names1))
clinical.data$stage <- mapvalues(clinical.data$stage, from = names1, to = names2)

names1 <- paste0("Stage II",c("","A", "B", "C"))
names2 <- rep("stage ii", length(names1))
clinical.data$stage <- mapvalues(clinical.data$stage, from = names1, to = names2)

names1 <- paste0("Stage III",c("","A", "B", "C"))
names2 <- rep("stage iii", length(names1))
clinical.data$stage <- mapvalues(clinical.data$stage, from = names1, to = names2)

names1 <- paste0("Stage IV",c("","A", "B", "C"))
names2 <- rep("stage iv", length(names1))
clinical.data$stage <- mapvalues(clinical.data$stage, from = names1, to = names2)

# all others categories
names1 <- c("[Not Applicable]", "[Not Available]", "[Discrepancy]", "[Unknown]", "Stage X") # not applicable is tissue type dependent
names2 <- rep("stage x", length(names1))
clinical.data$stage <- mapvalues(clinical.data$stage, from = names1, to = names2)

write.table(clinical.data, file="tcga-clinical-information.txt", quote=FALSE, sep="\t", row.names=FALSE)

