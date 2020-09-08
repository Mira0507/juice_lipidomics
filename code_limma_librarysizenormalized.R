library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
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

names(anno) <- str_replace_all(colnames(anno), " ", "_")


inner_join(tb, anno, by = c("Mass" = "mass_spectrum"))
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

############################################ Creating ExpressionSet Object ############################################ 

# Creating an ExpressionSet object
library(Biobase)
eset <- ExpressionSet(assayData = esiPN_mtx1,
                      phenoData = AnnotatedDataFrame(pheno),
                      featureData = AnnotatedDataFrame(feature))

# distribution of counts across the groups and metabolites
boxplot(exprs(eset)[999,] ~ pData(eset)[, "Group"],
        main = fData(eset)[999, "Metabolite"])

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
        ggplot(aes(x = Group, 
                   y = metabolite, 
                   fill = Group)) + 
        geom_boxplot() + 
        theme_bw() +
        ylab("Count") + 
        ggtitle("Library Disze Distribution") 



plotDensity_fn <- function(data, title) {
        
        plotDensities(data, 
                      group = pData(eset)[, "Group"],
                      legend = "topright",
                      main = title)
}

# Library Size normalization 
# Prevent missing values by adding 0.1 to every value (removing zero count)
exprs(eset) <- (exprs(eset) + 0.1) / colSums(exprs(eset))
plotDensity_fn(eset, "Distribution of Normalized Counts")

# Log transformation
exprs(eset) <- log(exprs(eset))
plotDensity_fn(eset, "Distribution of Log-transformed Counts")



# check whether there is any missing value 
sum(is.nan(rowMeans(exprs(eset))))
sum(is.na(rowMeans(exprs(eset))))
sum(is.infinite(rowMeans(exprs(eset))))


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


# Distribution of adjusted p-val
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


VolcanoPlot_fn <- function(dt, tit) {
        
        ggplot(dt, 
               aes(x = logFC,
                   y = LogOdds,
                   color = Change)) + 
                geom_point(alpha = 0.5) +
                theme_bw() + 
                ggtitle(paste("Volcano Plot:", tit)) +
                xlab("Log2 Fold Change") +
                ylab("-Log10 Adjusted P-value") + 
                scale_color_manual(values = c("red", "#999999", "blue")) + 
                xlim(-10, 10) 
        
}

LFC_NestedTable <- LFC_Table %>%
        
        # nesting
        nest(-Juice) %>%
        
        # calculating logodds
        mutate(data = map(data,
                          ~ mutate(.x, 
                                   LogOdds = -log10(adj.P.Val))),
               
               # Determining insignificant vs significant metabolites 
               data = map(data,
                          ~ mutate(.x, 
                                   Change = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Up-regulated",
                                                      adj.P.Val < 0.05 & logFC < 0 ~ "Down-regulated",
                                                      adj.P.Val > 0.05 ~ "Insignificant"))),
               
               data = map(data, ~ mutate(.x,
                                         Change = factor(Change, 
                                                         levels =  c("Up-regulated",
                                                                     "Insignificant",                                                       
                                                                     "Down-regulated"))))) %>%
        # Creating volcano plots
        mutate(volcano = map2(data,
                              Juice,
                              ~ VolcanoPlot_fn(.x, .y))) 

# volcano plots
VolcanoPlots <- grid.arrange(LFC_NestedTable$volcano[[1]],
                             LFC_NestedTable$volcano[[2]],
                             ncol = 1)

# Effect of apple or cranberry juice in metabolite change 
limma_test <- as.data.frame(summary(decideTests(fit2)))
names(limma_test) <- c("Change", "Contrast", "Metabolite_Number")

# plotting 
MetaboliteNumber_plot <- limma_test[1:6, ] %>%
        
        # data cleaning
        mutate(Contrast = str_replace_all(Contrast, "Effect_", ""),
               Change = case_when(Change == "Up" ~ "Up-regulated",
                                  Change == "Down" ~ "Down-regulated",
                                  Change == "NotSig" ~ "Insignificant")) %>% 
        
        # plotting
        ggplot(aes(x = Contrast,
                   y = Metabolite_Number,
                   fill = Change,
                   color = Change)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_label_repel(aes(label = Metabolite_Number,
                             color = Change),
                         fill = "white",
                         size = 5,
                         position = position_dodge(width = 1)) + 
        theme_bw() + 
        theme(axis.title.x = element_blank()) + 
        ylab("Metabolite Number") +
        ggtitle("Effect of Apple or Cranberry Juice in Metabolite Change") 
        
          




LFC_NestedSigTable <- LFC_NestedTable %>%
        mutate(data = map(data, ~ filter(.x, Change != "Insignificant")),
               data = map(data, ~ arrange(.x, desc(logFC))))

LFC_NestedSigTable$data[[1]]$KEGG_ID

###################################### Dimensionality Reduction & Clustering ###################################### 



################### PCA
library(factoextra)

# Running PCA
pca <- prcomp(exprs(eset),
              scale = TRUE, 
              center = TRUE)

# biplots
fviz_pca_var(pca, 
             col.var = "contrib",
             gradient.cols = c("blue", "red"),
             repel = TRUE, 
             title = "PCA")

fviz_pca_biplot(pca, 
                repel = TRUE,
                geom.ind = "point",
                title = "PCA")


# Computing cumulative proportion of var explained
v <- pca$sdev^2
pve <- v / sum(v)


pca_dt <- data.table(Prop_var = pve[1:10],
                 Principal_Component = 1:10)

# scree plot
PCA_Scree <- ggplot(pca_dt, aes(x = Principal_Component, 
                                y = Prop_var)) + 
        geom_line(size = 1, color = "blue") + 
        geom_point(size = 1.5, color = "blue") + 
        theme_bw() +
        scale_x_continuous(n.breaks = 10) + 
        ylab("Proportion of Variance Explained (%)") + 
        ggtitle("Scree Plot") + 
        xlab("Principal Component (PC)") +
        geom_vline(xintercept = 2, 
                   color = "red", 
                   size = 1)

# heatmap
pheatmap(pca$x[,1:2], 
         main = "PCA-heatmap")

# Hierarcical Clustering
PCA_hclustering <- cutree(hclust(dist(pca$x[,1:2]), 
                             method = "average"),
                      k = 3)

#Plotting
PCA_hClustering_Plot <- as.data.frame(pca$x[,1:2]) %>%
        rownames_to_column(var = "Metabolite") %>%
        mutate(hClustering = factor(PCA_hclustering)) %>%
        ggplot(aes(x = PC1, 
                   y = PC2, 
                   color = hClustering)) + 
        geom_point(alpha = 0.5, size = 2) + 
        theme_bw() + 
        ggtitle("PCA and Hierarchical Clustering")


################### t-SNE 


# Running tSNE
library(Rtsne)
tSNE_fn <- function(dt, pp) {
        
        # tSNE
        set.seed(2377)
        Rtsne(as.matrix(dt[, 2:ncol(dt)]),
              PCA = T, 
              perplexity = pp,
              max_iter = 2000)
}


tSNE_pp5 <- tSNE_fn(exprs(eset), 5)
tSNE_pp7 <- tSNE_fn(exprs(eset), 7)
tSNE_pp10 <- tSNE_fn(exprs(eset), 10)
tSNE_pp25 <- tSNE_fn(exprs(eset), 25)


plot(x = tSNE_pp5$Y[, 1],
     y = tSNE_pp5$Y[, 2])

plot(x = tSNE_pp7$Y[, 1],  # perplexity = 7 seems the best 
     y = tSNE_pp7$Y[, 2])

plot(x = tSNE_pp10$Y[, 1],
     y = tSNE_pp10$Y[, 2])

plot(x = tSNE_pp25$Y[, 1],
     y = tSNE_pp25$Y[, 2])

# heatmap
pheatmap(tSNE_pp7$Y, 
         main = "tSNE-heatmap")

# Hierarcical Clustering
tSNE_hclustering <- cutree(hclust(dist(tSNE_pp7$Y), 
                                 method = "average"),
                          k = 3)

tSNE_df <- as.data.frame(tSNE_pp7$Y)
rownames(tSNE_df) <- rownames(exprs(eset))
colnames(tSNE_df) <- c("Dim1", "Dim2")

# Plotting 
tSNE_hClustering_Plot <- tSNE_df %>%
        rownames_to_column(var = "Metabolite") %>%
        mutate(hClustering = factor(tSNE_hclustering)) %>%
        ggplot(aes(x = Dim1, 
                   y = Dim2, 
                   color = hClustering)) + 
        geom_point(alpha = 0.5, size = 2) + 
        theme_bw() + 
        ggtitle("tSNE and Hierarchical Clustering")