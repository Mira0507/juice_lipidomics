library(data.table)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(omu)

# Importing metabolite annotation 
anno <- rbind(fread("positive_lipid.txt"),
              fread("negative_lipid.txt")) 

names(anno) <- str_replace_all(colnames(anno), " ", "_")


# Creating count data frame
esiPN <- rbind(fread("ST000292_AN000466_ESI_POSITIVE_mwTab.txt"),
               fread("ST000292_AN000467_ESI_NEGATIVE_mwTAB.txt"))[-1, ] 

esiPN <- esiPN[!duplicated(esiPN$Samples),]

esiPN <- cbind(esiPN[, 1], 
               apply(esiPN[, -1], 2, as.numeric) + 0.1)

esiPN1 <- right_join(anno[, c("metabolite_name", "KEGG_ID")], 
                     esiPN, 
                     by = c("metabolite_name" = "Samples")) %>%
        filter(!is.na(KEGG_ID), KEGG_ID != "") %>% 
        rename(Metabolite = metabolite_name,
               KEGG = KEGG_ID)
esiPN2 <- esiPN1[!duplicated(esiPN1$Metabolite),]

# Creating meta data
Meta <- data.table(Sample = colnames(esiPN2[, -c(1, 2)])) %>%
        mutate(Group = case_when(str_detect(Sample, "a") ~ "Apple_Juice",
                                 str_detect(Sample, "b") ~ "Control",
                                 str_detect(Sample, "c") ~ "Cranberry_Juice"),
               Group = factor(Group, 
                              levels = c("Control",
                                         "Apple_Juice",
                                         "Cranberry_Juice")))

# Hierarchical class data 
DF <- assign_hierarchy(count_data = esiPN2,
                       keep_unknowns = TRUE,
                       identifier = "KEGG")

# Calculating fold change 
DF_Apple <- omu_summary(count_data = DF,
                      metadata = Meta,
                      response_variable = "Metabolite",
                      denominator = "Control", # control group
                      numerator = "Apple_Juice", # experimental group
                      Factor = "Group",
                      log_transform = T, 
                      p_adjust = "BH",
                      test_type = "welch")

DF_Cranberry <- omu_summary(count_data = DF,
                            metadata = Meta,
                            response_variable = "Metabolite",
                            denominator = "Control",
                            numerator = "Cranberry_Juice",
                            Factor = "Group",
                            log_transform = T, 
                            p_adjust = "BH",
                            test_type = "welch")


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

