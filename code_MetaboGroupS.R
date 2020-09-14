library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(limma)

# Importing metabolite annotation 
anno <- rbind(fread("positive_lipid.txt"),
              fread("negative_lipid.txt")) 

names(anno) <- str_replace_all(colnames(anno), " ", "_")



############################################ Data Cleaning ############################################ 


# Importing count matrix
esiPN <- rbind(fread("ST000292_AN000466_ESI_POSITIVE_mwTab.txt"),
               fread("ST000292_AN000467_ESI_NEGATIVE_mwTAB.txt"))[-1, ] %>%
        rename(metabolite_name = Samples)

# Converting to integer matrix
esiPN_mtx <- apply(as.matrix(esiPN[, -1]), 2, as.integer)

# naming the rows
rownames(esiPN_mtx) <- esiPN$metabolite_name

# removing duplicated rows and rows with NAs
esiPN_mtx1 <- esiPN_mtx[complete.cases(esiPN_mtx) &
                                !duplicated(rownames(esiPN_mtx)),]


# Creating phenotype data
pheno <- data.table(Sample = colnames(esiPN_mtx1)) %>%
        mutate(Group = case_when(str_detect(Sample, "a") ~ "Apple_Juice",
                                 str_detect(Sample, "b") ~ "Control",
                                 str_detect(Sample, "c") ~ "Cranberry_Juice"),
               Group = factor(Group, 
                              levels = c("Control",
                                         "Apple_Juice",
                                         "Cranberry_Juice"))) %>%
        column_to_rownames(var = "Sample") 

# rownames(Meta) and colnames(esiPN_mtx) have to be the same        
identical(rownames(pheno), colnames(esiPN_mtx1))

# Creating feature data 
anno1 <- anno[anno$metabolite_name %in% rownames(esiPN_mtx1),]
anno2 <- anno1[!duplicated(anno1$metabolite_name),]

feature <- data.frame(Metabolite = rownames(esiPN_mtx1)) %>%
        left_join(anno2, by = c("Metabolite" = "metabolite_name")) %>%
        column_to_rownames(var = "Metabolite")

# Creating Peaks Data 
write.csv(as.data.frame(esiPN_mtx1), file = "PeakData.csv")

# Creating sample data 
SampleData <- data.frame(sample = rownames(pheno),
                         batch = "1",
                         class = pheno$Group,
                         order = 1:51)
write.csv(SampleData, file = "SampleData.csv")
