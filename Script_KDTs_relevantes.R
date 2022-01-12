##############################################################################
################## Obtención KDTs relevantes y fármacos  #####################
######## Análisis de enriquecimiento funcional de esos KDTs relevantes #######
########################### ANA LAGO SAMPEDRO ################################
##############################################################################

# Establecer el directorio de trabajo (donde se guardarán los resultados)
setwd('~path/RESULTADOS')

#Cargar librerías necesarias
pacman::p_load("here", "rentrez","reutils","hipathia", "biomaRt", "utils", 
               "stringr", "AnnotationDbi", "org.Hs.eg.db",
               "dplyr","tidyr", "openxlsx", "data.table", "scales", 
               "BiocManager", "NMF", "clusterProfiler")

### modelos SHAP ########
data_folder = "~path/FILES/ML"
#Cargar los resultados de estabilidad para todos los circuitos procedentes del mapa
stability <- fread(file = file.path(data_folder, 
                                    "performance_stability_results_symbol.tsv")) %>% as.data.frame(.)
#View(stability)
rownames(stability) <-stability$V1
stability <- stability[-1,-1]
stability[is.na(stability)] = 0

#Cargar la matriz de puntuaciones de relevancias con la matriz delumbral de selección para filtrar
shap_entrez <-  fread(file = file.path(data_folder,"shap_summary_symbol.tsv"), 
                      header = T) %>% as.data.frame()
#View(shap_entrez)
rownames(shap_entrez)<- shap_entrez$V1
shap_entrez <- shap_entrez[ ,-1] 
shap_entrez_stable <- shap_entrez[rownames(shap_entrez) %in% rownames(stability)[stability$stability >= 0.4],]
#View(shap_entrez_stable)
threshold_entrez <- fread(file = file.path(data_folder,"shap_selection_symbol.tsv"), 
                          header = T) %>% as.data.frame() 
#View(threshold_entrez)
rownames(threshold_entrez) <- threshold_entrez$V1
threshold_entrez <- threshold_entrez[ ,-1]
threshold_entrez_stable <-  threshold_entrez[rownames(threshold_entrez) %in% rownames(stability)[stability$stability >= 0.4],]
#View(threshold_entrez_stable)

#Subset solo los valores shap que son relevantes para al menos un circuito en la matriz umbral
shap_relevant_stable <- shap_entrez_stable[,(apply(threshold_entrez_stable, 2, function(y) any(y == 1)))]
dim(shap_relevant_stable)
#View(shap_relevant_stable2)
shap_relevant_stable2<-cbind(rownames(shap_relevant_stable),shap_relevant_stable)
write.table(shap_relevant_stable2, "Shap_relevant_stable.tsv", row.names = FALSE, sep = "\t")

#Primero hay que reescalar por circuito (filas) cada matriz SHAP para hacer el circuito comparable
mat <- as.data.frame(t(apply(shap_relevant_stable, 1, function(x) x/max(abs(x))))) 
#Así, reescala los valores de puntuación para que estén en escala -1,1
#View(mat)

#Con esta matriz mat se hace el heatmap:
#con paquete NMF
aheatmap(mat, scale = "col", legend = TRUE, cellwidth = 10, cellheight = 5, 
         color = colorRampPalette(c("red","white","blue"))(50),
         Rowv = NA, fontsize = 8, filename = "aheatmap.png")
#unlink('aheatmap.pdf')


#Ahora con los KDTs de entrez o Symbols puedo buscar las drogas en el .tsv de DrugBank:
#Ordenarla para sacar la lista de KDTs que afectan de más circuitos a menos:
## Subset  solo los valores Shap que son relevantes para al menos un circuito.
KDTs <- threshold_entrez_stable[,(apply(threshold_entrez_stable, 2, function(y) any(y == 1)))]
dim(KDTs) #155 circuitos y 121 KDTs relevantes
#View(KDTs2)
KDTs2<-cbind(rownames(KDTs),KDTs)
write.table(KDTs2, "Circuits_y_KDTs.tsv", row.names = FALSE, sep = "\t")


#Suma de circuitos a los que afecta cada KDT:
circuitos_KDTs<-as.data.frame(colSums(KDTs))
KDTsn<-rownames(circuitos_KDTs)
#View(as.data.frame(KDTsn))
circuits_KDTs<-cbind(KDTsn,circuitos_KDTs)
#View(circuits_KDTs)
KDTs_order<-circuits_KDTs[order(circuitos_KDTs$`colSums(KDTs)`,decreasing = TRUE),]
#View(KDTs_order)

#Buscar las drogas dirigidos a los KDTs relevantes obtenidos:
### DRUGBANK ########
drug_bank = "~path/FILES"
## Cargar el fichero con todos los KDTs y sus drogas ### 
drugs <- fread(file = file.path(drug_bank, "drugbank_drug-bindings_v5.1.8.tsv")) %>% as.data.frame(.)
View(drugs)
#Filtrar KDTs para fármacos aprobados de humanos:

#Cargar librerías necesarias;
pacman::p_load("biomaRt", "org.Hs.eg.db", "magrittr", 
               "AnnotationDbi", "regexPipes","here")

#Seleccionar aquellos fármacos aprobados y que no contengan el tag "Withdrawn"
drugbank_approved <- drugs[-(regexPipes::grep(drugs$groups,
                                                        "withdrawn")),] %>% .[(regexPipes::grep(.$groups,
                                                                                                "^approved|approved,")),] %>% .[(regexPipes::grep(.$drug_binding,
                                                                                                                                                  "target_gene")),]                                                                                             
                                                                                                
#La tubería %>% pasa la salida del marco de datos que resulta de la función justo 
#antes de la tubería para ingresarla como el primer argumento de la función justo 
#después de la tubería.
#Y el punto hace referencia a la base de datos drugbank_alltar, pero la que sale 
#del pipe anterior.
#chequear los filtros para comprobar que ha ido bien
which(table(drugbank_approved$groups) > 0)
which(table(drugbank_approved$drug_binding) > 0 )
#Seleccionar dianas aprobadas con acción farmacológica conocida
drugbank_app_action <- drugbank_approved[(grep(drugbank_approved$pharmacological_action,"yes")),] %>% .[.$organism == "Humans",] %>% .[-which(is.na(.$entrez_id)),]
#View(drugbank_app_action)
#Chequear otra vez que los filtros han ido bien
table(drugbank_app_action$pharmacological_action)
which(table(drugbank_app_action$organism) > 0)
which(is.na (drugbank_app_action$entrez_id))
#Seleccionar solo genes humanos
entrez_a_targets <- as.character(unique(drugbank_app_action$entrez_id))
length(entrez_a_targets)
entrez_a_targets <- unique((mapIds(org.Hs.eg.db, keys=entrez_a_targets, keytype ="ENTREZID", column="SYMBOL", multiVals="first")))
#View(drugbank_app_action)
dim(drugbank_app_action)
length(unique(drugbank_app_action$drug_id))  #1461
length(unique(drugbank_app_action$entrez_id)) #714
approved_KDTs <- drugbank_app_action[ ,c(1,2,8,13,14)]
#View(approved_KDTs)
is.data.frame(approved_KDTs)
#View(approved_KDTs)

#Buscar la lista de KDTs relevantes en el fichero de approved_KDTs (approved_KDTs$gene_name):
#Lista de genes KDTs relevantes para el mapa DM2 ()
list_KDTs<-rownames(KDTs_order)
list_KDTs

Drugs_KDTs <- NULL
for(i in 1:length(list_KDTs)){
  Drugs_KDTs <- rbind(Drugs_KDTs, approved_KDTs[approved_KDTs$gene_name == list_KDTs[i],])
}
View(Drugs_KDTs)
length(unique(Drugs_KDTs$gene_name)) #120 KDTs
length(unique(Drugs_KDTs$drug_id))  #348 fármacos

write.table(Drugs_KDTs, "Drugs_KDTs.tsv",sep = "\t", row.names = FALSE)
##############################################################################


#Enriquecimiento funcional genes KDT relevantes 
############Enriquecimiento funcional#############
###############ClusterProfilerR###################
#A partir de la lista de KDTs relevantes:
length(list_KDTs)
head(list_KDTs)
KDTs<-list_KDTs
#Tengo la lista de 121 genes KDT con los que quiero hacer enriquecimiento funcional ORA:

### GO enrichment analysis ###
#GO comprende 3 ontologias ortogonales (MF molecular function; BP biological process;
#CC cellular component).
#Los análisis GO (groupGO, enrichGO y gseGO) para organismos en el objeto OrgDb, como Humano.
#Si tenemos datos GOannotation en data.frame con primera columna ID de genes
#y segunda columna con ID de GO, se puede usar las funciones enricher() y gseGO() para
#obtener test de sobrerrepresentación y analisis de enriquecimiento de un set de genes.

#Clasificación GO:
#En clusterProfiler, la función groupGO() está diseñada pra clasificación de genes basada en 
#distribución GO a un nivel específico.
library(clusterProfiler)
library(org.Hs.eg.db)

#Con ClusterProfiler voy a hacer clasificacion GO de esos efectores:

### GO over-representation analysis (ORA) ###
#La función enrichGO() para hacer el test de sobrerrepresentación de gene ontology
# Convert gene IDs to ENTREZID
# We will lose some genes here because not all IDs will be converted
ids<-bitr(KDTs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
head(ids)
dim(ids)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
dim(dedup_ids)

#Analisis de sobrerrepresentación de GO
go_enrich <- enrichGO(gene = dedup_ids$ENTREZID, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "all",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE, pool = FALSE)
#View(as.data.frame(go_enrich))
## Upset Plot
#Emphasizes the genes overlapping among different gene sets.
#library(enrichplot)
#upsetplot(go_enrich)
##Barplot
library(ggplot2)
barplot(go_enrich, split="ONTOLOGY",
        drop = TRUE, showCategory = 8, 
        title = "GO Biological Pathways",
        font.size = 6, order=T)+facet_grid(ONTOLOGY~.,scale="free")
##Dotplot
dotplot(go_enrich, showCategory = 10, title = "GO Biological Proccess", font.size=5, orderBy = "GeneRatio", 
        split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")

###Análisis WikiPathways
#WikiPathways es una base de datos continuamente en actualización por la comunidad de investigadores.
#Organismos que soporta; 
#get_wp_organisms()
#Para el análisis de sobrerrepresentación con este;
wiki <- enrichWP(dedup_ids$ENTREZID, organism = "Homo sapiens")
head(wiki,20)
dim(wiki)
##Barplot WikiPathways
barplot(wiki, showCategory = 10, order = T,
        title = "Enriched Pathways on WikiPathways", font.size = 8)
## Dotplot WikiPathways
dotplot(wiki, showCategory = 10, 
        title = "Enriched Pathways on WikiPathways",
        font.size = 7, orderBy = "GeneRatio")

###REACTOME
#ReactomePA está diseñada para análisis basado en rutas curadas de reactoma. Emplea modelo hipergeométrico para el análisis de sobrerrepresentación de genes. Hace corrección FDR.
#BiocManager::install("ReactomePA")
library(ReactomePA)
react <- enrichPathway(gene = dedup_ids$ENTREZID, 
                       pvalueCutoff = 0.05, readable = TRUE)
head(react)
dim(react)
##Barplot REACTOME
barplot(react, showCategory = 20, order = T,
        title = "Enriched Pathways on REACTOME", font.size = 6)
## Dotplot REACTOME
dotplot(react, showCategory = 15, 
        title = "Enriched Pathways on REACTOME",
        font.size = 6, orderBy = "GeneRatio")

### KEGG Orthology enrichment analysis ###
search_kegg_organism("hsa", by="kegg_code")
Human <- search_kegg_organism("Homo sapiens", by="scientific_name")
dim(Human)
head(Human)
#Analisis de enriquecimiento KEGG:
#KEGG pathway over-representation analysis
#necesito los identificadores de entrez
kk <- enrichKEGG(gene = dedup_ids$ENTREZID,
                 organism = "hsa", pvalueCutoff = 0.05)
#Bar Plot
barplot(kk, showCategory=10,label_format = 30, font.size = 7, 
        title = "Bar Plot KEGG paths", order=T)
#Dot plot
#Semejante a bar plot con la capacidad de encode otros score como dot size.
dotplot(kk, showCategory=12, font.size= 7, order=T) + ggtitle("Enriched pathways on KEGG")










