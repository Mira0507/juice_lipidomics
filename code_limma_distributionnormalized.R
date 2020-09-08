library(data.table)
library(tidyverse)
library(ggplot2)
library(Rtsne)
library(pheatmap)
library(gridExtra)
library(limma)
library(LipidMS)
library(LipidMSdata)



fas <- c("8:0", "10:0", "12:0", "14:0", "14:1", "15:0", "16:0", "16:1",
         "17:0", "18:0", "18:1", "18:2", "18:3", "18:4", "20:0", "20:1", "20:2",
         "20:3", "20:4", "20:5", "22:0", "22:1", "22:2", "22:3", "22:4", "22:5",
         "22:6", "24:0", "24:1", "26:0")
sph <- c("16:0", "16:1", "18:0", "18:1")
dbs <- createLipidDB(lipid = "all", chains = fas, chains2 = sph)

tb <- dbs[[1]]
for (i in 2:length(dbs)) {
        
        
        tb <- rbind(tb, dbs[[2]])
}

# Importing metabolite annotation 
anno <- rbind(fread("positive_lipid.txt"),
              fread("negative_lipid.txt")) 

############################################ Data Cleaning ############################################ 


# Importing count matrix
esiPN <- rbind(fread("ST000292_AN000466_ESI_POSITIVE.txt"),
               fread("ST000292_AN000467_ESI_NEGATIVE.txt"))

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


############################################ Creating ExpressionSet Object ############################################ 

# Creating an ExpressionSet object
library(Biobase)
eset <- ExpressionSet(assayData = esiPN_mtx1,
                      phenoData = AnnotatedDataFrame(pheno),
                      featureData = AnnotatedDataFrame(feature))

# distribution of counts across the groups and metabolites
boxplot(exprs(eset)[3370,] ~ pData(eset)[, "Group"],
        main = fData(eset)[3370, "Metabolite"])

############################################ Pre-processing ############################################ 

# Comparing library size: If library size distribution across the groups 
# is consistent, do librarysize normalization. Do distribution normalization
# (e.g. quantile normalization), otherwise. 
pheno1 <- rownames_to_column(pheno, 
                             var = "Sample") %>%
        mutate(Group = factor(Group, 
                              levels = c("Control",
                                         "Apple_Juice",
                                         "Cranberry_Juice")))

LibrarySize <- data.frame(metabolite = colSums(exprs(eset))) %>%
        rownames_to_column(var = "Sample") %>% 
        inner_join(pheno1, by = "Sample")

LibrarySize_Plot <- LibrarySize %>%
        ggplot(aes(x = Sample, 
                   y = metabolite, 
                   color = Group)) + 
        geom_point(size = 2) + 
        theme_bw() +
        theme(axis.text.x = element_blank()) + 
        ylab("Count") + 
        ggtitle("Library Disze Distribution") 



plotDensity_fn <- function(data, title) {
        
        plotDensities(data, 
                      group = pData(eset)[, "Group"],
                      legend = "topright",
                      main = title)
}

# Library Size normalization 
colsum <- apply(exprs(eset), 2, colSums)
# Log transformation
exprs(eset) <- log(exprs(eset))
plotDensity_fn(eset, "Distribution of Log-transformed Counts")

# Quantile normalization
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensity_fn(eset, "Distribution of Normalized Counts")



# Filtering out metabolites whose normalized counts are missing values 
keep1 <- !is.nan(rowMeans(exprs(eset)))
keep2 <- !is.infinite(rowMeans(exprs(eset)))
eset <- eset[keep1 & keep2, ]
plotDensity_fn(eset, "Counts after Filtering")


############################################ Sample Inspection ############################################

# Normalized count heatmap 
pheatmap(exprs(eset),
         annotation = select(pData(eset), Group),
         main = "Heatmap of Normalized Counts")

# Normalized count MDS plot
plotMDS(eset, 
        labels = pData(eset)[, "Group"],
        main = "Principal Component Analysis")

# Correlation heatmap
CorMat <- cor(exprs(eset))
pheatmap(CorMat,
         annotation = select(pData(eset), Group),
         main = "Correlation Heatmap")



############################################ Contrasts ############################################.

# Creating a design matrix
design <- model.matrix(~ 0 + Group, data = pData(eset))

colnames(design) <- levels(pData(eset)$Group)

ConMat <- makeContrasts(Effect_AppleJuice = Apple_Juice - Control,
                        Effect_CranberryJuice = Cranberry_Juice - Control,
                        Interaction = Cranberry_Juice - Apple_Juice,
                        levels = design)


################################ Fitting Linear Regression Model ################################ 

fit1 <- lmFit(eset, design)
fit2 <- contrasts.fit(fit1, contrasts = ConMat)
fit2 <- eBayes(fit2)


###################################### Extracting Results ###################################### 


extract_res <- function(data, ce, jc) {
        
        # Extracting LFC
        ResTable <- topTable(data, 
                 coef = ce, 
                 number = nrow(data)) 
        
        # Data Cleaning
        ResTable1 <- ResTable %>%
                rownames_to_column(var = "Metabolite") %>%
                mutate(Juice = jc)
        
        return(ResTable1)
}

res_AJ <- extract_res(fit2, "Effect_AppleJuice", "Apple")
res_CJ <- extract_res(fit2, "Effect_CranberryJuice", "Cranberry")
res_both <- topTable(fit2, NULL)



LFC_Table <- rbind(res_AJ, res_CJ) 

PVal_Dist_Plot <- ggplot(LFC_Table, 
                         aes(x = adj.P.Val,
                             fill = Juice,
                             color = Juice)) + 
        geom_density(alpha = 0.5) + 
        facet_grid(.~ Juice) +
        theme_bw() +
        xlab("Adjusted P-value") +
        ylab("Density") +
        ggtitle("Distribution of P-values")


LFC_NestedTable <- LFC_Table %>%
        nest(-Juice) %>%
        mutate(data = map(data,
                          ~ mutate(.x, 
                                   LogOdds = -log10(adj.P.Val))),
               data = map(data,
                          ~ mutate(.x, 
                                   Change = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
                                                      adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
                                                      adj.P.Val > 0.05 ~ "Insignificant"))),
               data = map(data, 
                          ~ mutate(.x, 
                                   Change = factor(Change,
                                                   levels = c("Up-regulated", 
                                                              "Down-regulated", 
                                                              "Insignificant")))))


limma_test <- as.data.frame(summary(decideTests(fit2)))
names(limma_test) <- c("Change", "Contrast", "Metabolite_Number")

