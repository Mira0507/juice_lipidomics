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

LFC_NestedSigTable$data[[1]]


# Creating a heatmap for differentially abundant metabolites 
SigMetabolites <- rbind(exprs(eset)[LFC_NestedSigTable$data[[1]]$Metabolite, ],
                        exprs(eset)[LFC_NestedSigTable$data[[2]]$Metabolite, ])

SigMetabolites <- SigMetabolites[!duplicated(SigMetabolites), ]


pheatmap(SigMetabolites,
         annotation = select(pData(eset), Group),
         main = "Differentially Abundant Metabolites")
