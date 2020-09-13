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
        ggtitle("Library Size Distribution") 



plotDensity_fn <- function(data, leg, title) {
        
        plotDensities(data, 
                      group = pData(eset)[, "Group"],
                      legend = leg,
                      main = title)
}

# count distribution before normalization
plotDensity_fn(eset, 
               "topright",
               "Distribution of Counts before Log-transformation")

MDPlot_before_Preprocessing <- plotMD(eset,
                                      main = "MD Plot before Preprocessing",
                                      xlab = "Mean Log-count",
                                      ylab = "Log Fold Change")

# Library Size normalization 
# Prevent missing values by adding 0.1 to every value (removing zero count)
exprs(eset) <- (exprs(eset) + 1) / colSums(exprs(eset))
CounDist_after_Norm <- plotDensity_fn(eset, 
                                       "topright",
                                       "Distribution of Normalized Counts")

# Log transformation
exprs(eset) <- log(exprs(eset))
CounDist_after_LogTrans <- plotDensity_fn(eset, 
                                          "topleft",
                                          "Distribution of Log-transformed Counts")
# filtering out low count metabolites
keep <- rowMeans(exprs(eset)) > -15

# updating an ExpressionSet object
eset <- ExpressionSet(assayData = exprs(eset)[keep, ],
                      phenoData = AnnotatedDataFrame(pheno),
                      featureData = AnnotatedDataFrame(feature[keep, ]))
        
CounDist_after_Filt <- plotDensity_fn(eset, 
                                      "topleft",
                                      "Distribution of Filtered Counts")


# check whether there is any missing value 
sum(is.nan(rowMeans(exprs(eset))))
sum(is.na(rowMeans(exprs(eset))))
sum(is.infinite(rowMeans(exprs(eset))))


############################################ Sample Inspection ############################################

# MD plot
MDPlot_after_Preprocessing <- plotMD(eset,
                                     main = "MD Plot after Preprocessing",
                                     xlab = "Mean Log-count",
                                     ylab = "Log Fold Change")

# Normalized count heatmap 
CountHeatmap <- pheatmap(exprs(eset),
                         annotation = select(pData(eset), Group),
                         main = "Heatmap of Normalized Counts")

# Normalized count MDS plot
CountMDS <- plotMDS(eset, 
                    labels = pData(eset)[, "Group"],
                    col = as.numeric(pData(eset)[, "Group"]),
                    main = "Multidimensional Scaling (MDS) plot")

# Correlation heatmap
CorMat <- cor(exprs(eset))
CorHeatmap <- pheatmap(CorMat,
                       annotation = select(pData(eset), Group),
                       main = "Correlation Heatmap")


####################### Omu 
library(omu)

# Creating count data frame
Omu_esiPN <- rbind(fread("ST000292_AN000466_ESI_POSITIVE_mwTab.txt"),
                   fread("ST000292_AN000467_ESI_NEGATIVE_mwTAB.txt"))[-1, ] 

Omu_esiPN <- Omu_esiPN[!duplicated(Omu_esiPN$Samples),]

Omu_esiPN <- cbind(Omu_esiPN[, 1], 
                   apply(Omu_esiPN[, -1], 2, as.numeric) + 0.1)


Omu_esiPN1 <- right_join(anno[, c("metabolite_name", "KEGG_ID")], 
                         Omu_esiPN, 
                         by = c("metabolite_name" = "Samples")) %>%
        filter(!is.na(KEGG_ID), KEGG_ID != "") %>% 
        rename(Metabolite = metabolite_name,
               KEGG = KEGG_ID)
Omu_esiPN2 <- Omu_esiPN1[!duplicated(Omu_esiPN1$Metabolite),]

# Creating meta data
Meta <- data.table(Sample = colnames(Omu_esiPN2[, -c(1, 2)])) %>%
        mutate(Group = case_when(str_detect(Sample, "a") ~ "Apple_Juice",
                                 str_detect(Sample, "b") ~ "Control",
                                 str_detect(Sample, "c") ~ "Cranberry_Juice"),
               Group = factor(Group, 
                              levels = c("Control",
                                         "Apple_Juice",
                                         "Cranberry_Juice")))

# Hierarchical class data 
DF <- assign_hierarchy(count_data = Omu_esiPN2,
                       keep_unknowns = TRUE,
                       identifier = "KEGG")

# plots
Omu_PCA <- 
        PCA_plot(count_data = DF,
                 metadata = Meta,
                 variable = "Group",
                 color = "Group", 
                 response_variable = "Metabolite") +
        theme_bw() +
        ggtitle("PCA") 

Omu_boxplot <- plot_boxplot(count_data = DF,
                            metadata = Meta,
                            Factor = "Group",
                            aggregate_by = "Class",
                            log_transform = T,
                            response_variable = "Metabolite",
                            fill_list = c("grey", "pink","lightblue")) + 
        ggtitle("Metabolite Counts") + 
        theme(axis.text.x = element_blank())

Omu_heatmap <- plot_heatmap(count_data = DF,
                            metadata = Meta,
                            Factor = "Group",
                            aggregate_by = "Class",
                            log_transform = T,
                            response_variable = "Metabolite",
                            high_color = "red",
                            low_color = "blue") +
        ggtitle("Abundance Heatmap of Metabolites") +
        ylab("Metabolite Class") + 
        theme(axis.text.x = element_blank())

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
        geom_density(alpha = 0.7) + 
        facet_grid(.~ Juice) +
        theme_bw() +
        xlab("Adjusted P-value") +
        ylab("Density") +
        ggtitle("Distribution of P-values (False Discovery Rate)")


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
limma_test1 <- decideTests(fit2)
limma_test2 <- as.data.frame(summary(limma_test1))
names(limma_test2) <- c("Change", "Contrast", "Metabolite_Number")

# plotting 
set.seed(16)
MetaboliteNumber_plot <- limma_test2[1:6, ] %>%
        
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



venndiagram <- vennDiagram(limma_test1,
                           main = "Number of Significantly Changed Metabolites")


LFC_NestedSigTable <- LFC_NestedTable %>%
        mutate(data = map(data, ~ filter(.x, Change != "Insignificant")),
               data = map(data, ~ arrange(.x, desc(logFC))))

LFC_NestedSigTable$data[[1]]


# Creating a heatmap for differentially abundant metabolites 
SigMetabolites <- rbind(exprs(eset)[LFC_NestedSigTable$data[[1]]$Metabolite, ],
                        exprs(eset)[LFC_NestedSigTable$data[[2]]$Metabolite, ])

SigMetabolites <- SigMetabolites[!duplicated(SigMetabolites), ]


SigMetabolite_Heatmap <- pheatmap(SigMetabolites,
                                  annotation = select(pData(eset), Group),
                                  main = "Differentially Abundant Metabolites")

###################################### Dimensionality Reduction & Clustering ###################################### 
library(ggfortify)
library(cluster)

################### PCA
library(factoextra)

# Running PCA
pca <- prcomp(SigMetabolites,
              scale = TRUE, 
              center = TRUE)

# biplots
fviz_pca_var(pca, 
             col.var = "contrib",
             gradient.cols = c("blue", "red"),
             repel = TRUE, 
             title = "PCA")

Biplot <- fviz_pca_biplot(pca, 
                          repel = TRUE,
                          geom.ind = "point",
                          title = "Biplot")


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
PCA_Heatmap <- pheatmap(pca$x[,1:2], 
                        main = "PCA-heatmap")



# Data Cleaning
set.seed(490)
PCA_kmclustering_dt <- as.data.frame(pca$x[,1:2]) %>%
        rownames_to_column(var = "Metabolite") %>%
        mutate(PCA_kmClustering = factor(kmeans(pca$x[,1:2], 
                                                centers = 2)$cluster)) 



#Plotting
set.seed(490)
PCA_kmClustering_Plot <- 
        autoplot(stats::kmeans(PCA_kmclustering_dt[, 2:3], 2),
                 data = PCA_kmclustering_dt,
                 frame = TRUE) +
        theme_bw() + 
        ggtitle("PCA and K-Means Clustering") 

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


tSNE_pp5 <- tSNE_fn(SigMetabolites, 5)
tSNE_pp7 <- tSNE_fn(SigMetabolites, 7)
tSNE_pp10 <- tSNE_fn(SigMetabolites, 10)



plot(x = tSNE_pp5$Y[, 1],
     y = tSNE_pp5$Y[, 2])

plot(x = tSNE_pp7$Y[, 1],  # perplexity = 7 seems the best 
     y = tSNE_pp7$Y[, 2])

plot(x = tSNE_pp10$Y[, 1],  
     y = tSNE_pp10$Y[, 2])


# heatmap
tsne7DT <- tSNE_pp7$Y
colnames(tsne7DT) <- c("Dim_1", "Dim_2")
rownames(tsne7DT) <- rownames(SigMetabolites)
tSNE_Heatmap <- pheatmap(tsne7DT, 
                         main = "tSNE-heatmap")

# Data Cleaning
tSNE_df <- as.data.frame(tSNE_pp7$Y)
rownames(tSNE_df) <- rownames(SigMetabolites)
colnames(tSNE_df) <- c("Dim1", "Dim2")

# Data Cleaning
set.seed(490)
tSNE_kmClustering_dt <- tSNE_df %>%
        rownames_to_column(var = "Metabolite") %>%
        mutate(tSNE_kmclustering = factor(kmeans(tSNE_df, 
                                                 centers = 2)$cluster))

# Plotting (w ggfortify and cluster package)
set.seed(490)
tSNE_kmClustering_Plot <- 
        autoplot(stats::kmeans(tSNE_kmClustering_dt[, 2:3], 2),
                 data = tSNE_kmClustering_dt,
                 frame = TRUE) +
        theme_bw() + 
        ggtitle("tSNE and K-Means Clustering") +
        xlab("Dim1 (89.09%)") +
        ylab("Dim2 (10.91%)")


# Data Cleaning for comparing clustering between PCA and tSNE
PCA_tSNE_Table <- inner_join(PCA_kmclustering_dt,
                             tSNE_kmClustering_dt,
                             by = "Metabolite") 

CombinedCluster <- as.data.table(table(PCA_tSNE_Table$PCA_kmClustering, 
                                       PCA_tSNE_Table$tSNE_kmclustering))

colnames(CombinedCluster) <- c("PCA",
                               "tSNE",
                               "Metabolite_Number")

# Plotting for comparing clustering between PCA and tSNE
CombinedCluster_plot <- 
        ggplot(CombinedCluster,
               aes(x = PCA, y = tSNE, label = Metabolite_Number)) + 
        geom_tile(data = CombinedCluster, 
                  aes(fill = Metabolite_Number),
                  alpha = 0.7) + 
        geom_text(size = 15) +
        theme_bw() + 
        ggtitle("Number of Metabolites Clustered by PCA and tSNE") + 
        xlab("Cluster by PCA") +
        ylab("Cluster by tSNE")

