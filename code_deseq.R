library(data.table)
library(tidyverse)
library(ggplot2)
library(Rtsne)
library(pheatmap)
library(gridExtra)
library(DESeq2)
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


# Creating metadata
Meta <- data.table(Sample = colnames(esiPN_mtx1)) %>%
        mutate(Group = case_when(str_detect(Sample, "a") ~ "Apple_Juice",
                                 str_detect(Sample, "b") ~ "Control",
                                 str_detect(Sample, "c") ~ "Cranberry_Juice")) %>%
        column_to_rownames(var = "Sample") 

# rownames(Meta) and colnames(esiPN_mtx) have to be the same        
identical(rownames(Meta), colnames(esiPN_mtx1))


# Creating DESeq object 
des <- DESeqDataSetFromMatrix(countData = esiPN_mtx1,
                              colData = Meta,
                              design = ~ Group)
# Calculating size factors
des <- estimateSizeFactors(des)
sizeFactors(des)

# normalization
norm_counts <- counts(des, normalized = TRUE)

# log transformation
vsd <- vst(des, blind = TRUE)

# Correlation heatmap
vsd_mat <- assay(vsd)
vsd_cor <- cor(vsd_mat)

CorHeatmap <- pheatmap(vsd_cor, 
                       main = "Correlation Heatmap",
                       annotation = select(Meta, Group))

# PCA 
CountPCA <- plotPCA(vsd, intgroup = "Group") + 
        ggtitle("PCA (Normalized & Log-transformed Counts)") + 
        theme_bw()


# Running analysis
des <- DESeq(des)

# Dispersion model
dispersion <- plotDispEsts(des)

# Extracting Log fold change (LFC)
LFCTable_fn <- function(res) {
        
        # converting to data frame 
        res <- as.data.frame(res)
        
        # adding metabolite variable 
        res$Matabolite <- rownames(res)
        
        raw_result <- res %>% mutate(Significance = ifelse(padj <= 0.05, "O", "X"),
                       LogOdd = -log(padj)) 
        
        arranged_result <- raw_result %>% 
                filter(Significance == "O") %>% 
                arrange(desc(log2FoldChange))
        
        return(list(raw_result, arranged_result))
                
}


control_v_apple <- results(des, 
                           alpha = 0.5,
                           contrast = c("Group",
                                        "Apple_Juice",
                                        "Control"))

Vapple_LFC <- LFCTable_fn(control_v_apple)[[1]]

Vapple_LFC_arranged <- LFCTable_fn(control_v_apple)[[2]]

control_v_cranberry <- results(des, 
                           alpha = 0.5,
                           contrast = c("Group",
                                        "Cranberry_Juice",
                                        "Control"))

Vcranberry_LFC <- LFCTable_fn(control_v_cranberry)[[1]]

Vcranberry_LFC_arranged <- LFCTable_fn(control_v_cranberry)[[2]]


# volcano plots
Vapple_LFC_plot <- 
        ggplot(Vapple_LFC, 
               aes(x = log2FoldChange,
                   y = LogOdd,
                   color = Significance)) + 
        geom_point(alpha = 0.3) +
        theme_bw() + 
        ggtitle("Lipid Change: Control vs Apple Juice")

Vcranberry_LFC_plot <- 
        ggplot(Vcranberry_LFC, 
               aes(x = log2FoldChange,
                   y = LogOdd,
                   color = Significance)) + 
        geom_point(alpha = 0.3) +
        theme_bw() + 
        ggtitle("Lipid Change: Control vs Cranberry Juice")




############################################### t-SNE ############################################### 


# t-SNE analysis
tsne_fn <- function(dt, pp) {
        
        # tSNE
        set.seed(161)
        Rtsne(dt,
              PCA = T,
              perplexity = pp,
              max_iter = 2000)
        
}

# Testing perplexity = 5, 10, 15, 20
tsne5 <- tsne_fn(norm_counts, 5)
tsne10 <- tsne_fn(norm_counts, 10)
tsne15 <- tsne_fn(norm_counts, 15)
tsne20 <- tsne_fn(norm_counts, 20)


plot(tsne5$itercosts, 
     type = "l",
     ylab = "Total K-L Divergence Cost",
     xlab = "Gradient Descent (50 Steps Each)",
     main = "Optimal Number of Iterations")

plot(tsne15$itercosts, 
     type = "l",
     ylab = "Total K-L Divergence Cost",
     xlab = "Gradient Descent (50 Steps Each)",
     main = "Optimal Number of Iterations")

# Data cleaning for t-SNE plotting
tSNE_table <- cbind(tsne5$Y, 
                    tsne10$Y,
                    tsne15$Y,
                    tsne20$Y) %>%
        as.data.table()

names(tSNE_table) <- c("tsne5_X", "tsne5_Y",
                       "tsne10_X", "tsne10_Y",
                       "tsne15_X", "tsne15_Y",
                       "tsne20_X", "tsne20_Y")


# Heatmap with eigen values 
pheatmap(tSNE_table[, 1:2]) # perplexity = 5
pheatmap(tSNE_table[, 3:4]) # perplexity = 10
pheatmap(tSNE_table[, 5:6]) # perplexity = 15 (looks like the best) 
pheatmap(tSNE_table[, 7:8]) # perplexity = 20

# Data cleaning
tSNE_table1 <- tSNE_table %>%
        gather(X, Eigenvalue_X, ends_with("X")) %>%
        gather(Y, Eigenvalue_Y, ends_with("Y")) %>%
        mutate(X = str_replace_all(X, "_X", "")) %>%
        mutate(Y = str_replace_all(Y, "_Y", "")) %>%
        mutate(X = factor(X, levels = paste0("tsne", c(5, 10, 15, 20))),
               Y = factor(Y, levels = paste0("tsne", c(5, 10, 15, 20))))
        

tSNE_table2 <- tSNE_table1[tSNE_table1$X == tSNE_table1$Y, ]

# t-SNE visualization 
ggplot(tSNE_table2, 
       aes(x = Eigenvalue_X, 
           y = Eigenvalue_Y)) + 
        geom_point(alpha = 0.3) + 
        facet_grid(X ~.) +
        theme_bw() + 
        xlab("Dim 1") + 
        ylab("Dim 2") + 
        ggtitle("t-SNE Visualization with Perplexity Change")

# Determining k: scree plot
Scree_fn <- function(dt, tit, intercept) {
        
        # measuring total within ss 
        set.seed(1987)
        ttWithinss <- map_dbl(1:10, 
                              function(k) {
                                      km <- kmeans(x = dt, 
                                                   centers = k,
                                                   nstart = 25)
                                      km$tot.withinss})
        screeDT <- data.table(
                k = 1:10, 
                Total_Within_SS = ttWithinss
        )
        
        # creating a scree plot
        ggplot(screeDT,
               aes(x = k,
                   y = Total_Within_SS)) + 
                geom_line(size = 1, color = "blue") +
                geom_point(size = 2, color = "blue") +
                ggtitle(tit) +
                theme_bw() + 
                ylab("Total Within Cluster Sum of Squares") + 
                scale_x_continuous(breaks = 1:10, 
                                   minor_breaks = NULL) +
                geom_vline(xintercept = intercept, 
                           size = 1,
                           color = "red")
}

tSNE_table3 <- tSNE_table2 %>%
        nest(Eigenvalue_X, Eigenvalue_Y) %>%  
        mutate(scree = map(data, ~ Scree_fn(.x, "Scree Plot", 4)),
               Distance = map(data, ~ dist(.x)),
               Clustering = map(Distance, ~ hclust(.x, method = "average")),
               Cluster = map(Clustering, ~ factor(cutree(.x, k = 4))),
               Clustering_Plot = map2(data, 
                                      Cluster, 
                                      ~ ggplot(.x, 
                                               aes(x = Eigenvalue_X,
                                                   y = Eigenvalue_Y,
                                                   color = .y)) + 
                                              geom_point(alpha = 0.2) + 
                                              theme_bw() +
                                              xlab("Dim 1") + 
                                              ylab("Dim 2")),
               HeatMap = map(data, ~ pheatmap(.x)))

grid.arrange(tSNE_table3$Clustering_Plot[[1]] + 
                     ggtitle("Hierarchical Clustering (Perplexity = 5)"),
             tSNE_table3$Clustering_Plot[[2]] + 
                     ggtitle("Hierarchical Clustering (Perplexity = 10)"),
             tSNE_table3$Clustering_Plot[[3]] + 
                     ggtitle("Hierarchical Clustering (Perplexity = 15)"),
             tSNE_table3$Clustering_Plot[[4]] + 
                     ggtitle("Hierarchical Clustering (Perplexity = 20)"),
             ncol = 2, nrow = 2)


