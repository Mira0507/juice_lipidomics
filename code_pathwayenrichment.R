###################################### Pathway Enrichment Analysis ###################################### 

#################### FELLA


library(FELLA)



PathwayEnrichment_fn <- function(keggID, 
                                 directory) {
        
        # building a graph 
        set.seed(11)
        gr <- buildGraphFromKEGGREST(organism = "hsa", 
                                     filter.path = keggID)
        
        # assigning directory 
        tempdir <- paste0(tempdir(), directory)
        unlink(tempdir, recursive = TRUE)
        
        
        buildDataFromGraph(keggdata.graph = gr,
                           databaseDir = tempdir, 
                           internalDir = FALSE,
                           matrices = "diffusion",
                           normality = "diffusion",
                           niter = 50)
        
        FD <- loadKEGGdata(databaseDir = tempdir,
                           internalDir = FALSE,
                           loadMatrix = "diffusion")
        
        return(FD)
        
        
        
}

# Fella: temporary local database 
FellaData_Apple <- PathwayEnrichment_fn(LFC_NestedSigTable$data[[1]]$KEGG_ID, 
                                        "/my_database1")

FellaData_Cranberry <- PathwayEnrichment_fn(LFC_NestedSigTable$data[[2]]$KEGG_ID, 
                                            "/my_database2")

# Creating an FELLA.USER object
PathwayObject_fn <- function(cpd, fella) {
        
        # Creating an analysis object
        analysis <- defineCompounds(cpd,
                                    data = fella)
        
        # approx = how to process statistics 
        # normality (normalization), simulation (empirical p-value), 
        # t (t-distribution), gamma (gamma distribution)
        analysis <- runDiffusion(object = analysis,
                                 data = fella,
                                 approx = "normality") 
        
        return(analysis) 
}

# check p-scores under 0.05  
Analysis_Ap <- PathwayObject_fn(LFC_NestedSigTable$data[[1]]$KEGG_ID,
                                FellaData_Apple)

Anakysis_Cr <- PathwayObject_fn(LFC_NestedSigTable$data[[2]]$KEGG_ID,
                             FellaData_Cranberry)
# mapped Analysis_Ap 
getInput(Analysis_Ap)
getInput(Anakysis_Cr)
# unmapped identifier 
getExcluded(Analysis_Ap)
getExcluded(Anakysis_Cr)

# Pathway map
jpeg("FellaMap_Apple.jpeg",
     width = 500, 
     height = 500, quality = 1000)
plot(Analysis_Ap, 
     method = "diffusion",
     data = FellaData_Apple,
     nlimit = 150,
     vertex.label.cex = 0.5)
dev.off()


PathwayTable_fn <- function(cpd, fella, analysis) {
        
        
        # nlimit and vertex.label.cex are tunable
        path_table <- generateResultsTable(object = analysis,
                                           method = "diffusion",
                                           nlimit = 150,
                                           data = fella)
        
        enzyme_table <- generateEnzymesTable(object = analysis,
                                             method = "diffusion",
                                             nlimit = 150,
                                             data = fella)
        
        return(list(path_table, enzyme_table))
        
}



PathwayRes_Apple <- PathwayTable_fn(LFC_NestedSigTable$data[[1]]$KEGG_ID,
                                    FellaData_Apple,
                                    Analysis_Ap)

PathwayRes_Cranberry <- PathwayTable_fn(LFC_NestedSigTable$data[[2]]$KEGG_ID,
                                        FellaData_Cranberry,
                                        Anakysis_Cr)


PathwayDataCleaning_fn <- function(dt) {
        
        data <- dt %>% filter(Entry.type == "pathway") %>% 
                arrange(p.score) %>%
                select(- Entry.type)
        
        names(data) <- c("KEGG_PathwayID", 
                         "Pathway_Name",
                         "P_Value")
        
        return(data)
                
}

KEGGPath_Apple <- PathwayDataCleaning_fn(PathwayRes_Apple[[1]])
KEGGPath_Cranberry <- PathwayDataCleaning_fn(PathwayRes_Cranberry[[1]])

PathTable_Apple <- formattable(KEGGPath_Apple)
formattable(KEGGPath_Cranberry) # There is no pathway significantly found (p <= 0.05)

#################### pathview


library(pathview)

# import required databases (from the pathview package)
data(cpd.accs)
data(cpd.names)
data(kegg.met)
data(ko.ids)
data(rn.list)
data(gene.idtype.list)
data(gene.idtype.bods)
data(cpd.simtypes)

LFC_NestedSigTable$data[[1]]
LFC_NestedSigTable$data[[2]]




# Creating a FC vector named by KEGG compound ID (descending order)

# C03340     C00989     C00741     C12284     C00109   
# 6.3300292  0.9735598  0.4220332  0.3251592  0.2406034 

AppleLFC <- LFC_NestedSigTable$data[[1]]$logFC
names(AppleLFC) <- LFC_NestedSigTable$data[[1]]$KEGG_ID
AppleLFC <- sort(AppleLFC, decreasing = T)


# KEGG pathview
# pathway.id = enriched pathways suggested by FELLA
pathview(cpd.data = AppleLFC,
         species = "hsa",
         pathway.id = KEGGPath_Apple$KEGG_PathwayID)


