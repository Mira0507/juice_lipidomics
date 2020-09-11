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