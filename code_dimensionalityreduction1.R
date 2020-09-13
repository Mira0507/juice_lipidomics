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

# Hierarcical Clustering
PCA_hclustering <- cutree(hclust(dist(pca$x[,1:2]), 
                                 method = "average"),
                          k = 2)


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

# Hierarcical Clustering
tSNE_hclustering <- cutree(hclust(dist(tsne7DT), 
                                  method = "average"),
                           k = 2)

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
        xlab("Dim1 (89.1%)") +
        ylab("Dim2 (10.9%)")



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
           
