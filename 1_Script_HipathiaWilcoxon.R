##############################################################################
############# ANALISIS DE DATOS DE ACTIVACIÓN CIRCUITOS RNASeq ###############
#### PARA ANALISIS DIFERENCIAL DEL PROYECTO GSE164416 (CON TEST WILCOXON) ####
########################### ANA LAGO SAMPEDRO ################################
##############################################################################

# Establecer el directorio de trabajo (donde se guardarán los resultados)
setwd('~path/RESULTADOS')

## ------------------------------------------------------------------------
# Instalar los paquetes necesarios:
BiocManager::install("hipathia")
BiocManager::install("edgeR")

#Cargar librerías:
library(edgeR)
library(hipathia)

## ------------------------------------------------------------------------
#Input: datos crudos del proyecto RNA-Seq en cuentas totales para normalizar
#y Metadatos (diseño del estudio indicando cada muestra y su condición)
## Matriz de counts:  counts_matrix(T2DvsND).txt
## Matriz de diseño:  design_matrix(T2DvsND).txt

## ------------------------------------------------------------------------
###### Definir directorios; matriz de conteos y matriz de diseño ######
data_folder<-"~path/FILES"

fname<-"counts_matrix(T2DvsND).txt"
fpath<-file.path(data_folder,fname)

dname<-"design_matrix(T2DvsND).txt"
dpath<-file.path(data_folder,dname)

###### Leer los datos de la matriz de counts
expreset_raw<-read.table(fpath,header=TRUE,sep="\t")
dim(expreset_raw)  
rownames(expreset_raw)<-expreset_raw$ensembl  #Asignar a cada fila el nombre del gen
expreset<-expreset_raw[,-1]  #Eliminar la primera columna, para tener una matriz 
#de counts exclusivamente
hhead(expreset)

###### Leer datos de Diseño del experimento
matrix_design <- read.table(dpath,header=FALSE,sep="\t")
names (matrix_design) = c("ID", "Condition","name")
rownames(matrix_design)<-matrix_design$ID  #Asignar a cada fila el nombre de la muestra
matrix_design<-matrix_design[,-1]  #Eliminar la primera columna
hhead(matrix_design)

#######################################################################################

###### Construir objeto DGE de EdgeR para la normalización de los datos con esta 
# librería y añadir diseño al DGE.
#counts= matriz de conteos
#group= información sobre el grupo/condición experimental para cada muestra 
# (T2D vs ND == Condition)
dge <- DGEList(counts = expreset, group = matrix_design$Condition)
hhead(dge)

###### Normalizar datos crudos por método TMM con calcNormFactors
tmm<-calcNormFactors(dge,method="TMM") #calcula el factor de normalización x TMM de edgeR
hhead(tmm)

###### Aplicar la normalización; añadir 3 counts a cada dato para después corregir 
# por el factor de normalización y a ese resultado se le hace transformación log2
logcpm<-cpm(tmm,prior.count=3,log=TRUE)
hhead(logcpm)

######################################################################################

###### HIPATHIA: Obtención de valores de activación de circuitos y posteriormente,
# se realiza el análisis de activación diferencial de circuitos entre DM2 vs ND.
# Se usa la información de las rutas fisiológicos de humano, que incluye 79 rutas,
# sin rutas de enfermedades.

# Primero, traducir los identificadores de genes de la matriz normalizada al 
 #específico para Hipathia (obtiene los IDs de Entrez)
trans_data<-translate_data(logcpm,"hsa")
###Resultado en el proyecto GSE164416###
# translated ids = 23982 (0.41) 
# untranslated ids = 34354 (0.59) 
# multihit ids = 185 (0.0032) 

# Escalar los valores al rango [0,1] para trabajar en Hipathia
exp_data<-normalize_data(trans_data)
head(exp_data)
dim(exp_data) #Se obtienen 23881 filas (genes) y 57 columnas (muestras)

### Calcular los valores de activación de circuitos en 79 rutas fisiológicas
# Cargar las rutas fisiológicas
physiological_path_list<-read.table("~path/physiological_paths.tsv",
                                    sep = "\t", stringsAsFactors = F,
                                    header = F, quote = "")
metaginfo<-load_pathways("hsa",pathways_list = physiological_path_list[,2])
### Obtener información del objeto metaginfo de hipathia conseguido con la 
# función load_pathways:
dim(as.data.frame(metaginfo$all.genes))  #incluye 2516 genes
dim(as.data.frame(metaginfo$path.norm))  #6546 circuitos con efectores y 79 rutas

# Calcular los valores de activación de rutas
results_physiologicals<-hipathia(exp_data,metaginfo,decompose=FALSE,verbose=FALSE)
###Resultado###
# Added missing genes: 24 (0.1%)

#La siguiente función devuelve el objeto con los niveles de activación de cada 
#para cada muestra:
node_data<-get_nodes_data(results_physiologicals, matrix = TRUE)
head(node_data)
dim(node_data)
#3944 nodos que aparecen en hipathia

#La siguiente función Obtiene la matriz de valores de activación de circuitos 
#para cada muestra:
path_vals<-get_paths_data(results_physiologicals,matrix=TRUE)
head(path_vals)
dim(path_vals)
#1098 circuitos que aparecen en hipathia


### Análisis de activación diferencial de circuitos 

# Primero se define el factor condición;
sample_group <- matrix_design[colnames(path_vals),"Condition"]
sample_group

# Hacer test estadístico de Wilcoxon
comp <- do_wilcoxon(path_vals, sample_group, g1="T2D", g2="ND", 
                    paired = FALSE, adjust = TRUE)
#el g2 es el grupo de referencia  (,order=TRUE pone los resultados por orden p.adjust)
head(comp)
top_pathways(comp)
#Guardar el data.frame comp para usar después en el fichero circuits.rds:
comp2<-cbind(rownames(comp),comp)
head(comp2)
# ¿Cómo obtener el nombre del circuito-efector? 
path_names <- get_path_names(metaginfo, rownames(comp))
head(path_names)
#Añadir los nombres de los circuitos
comp <- cbind(path_names, comp)
head(comp)
dim(comp)
#### ESCRIBIR ESTA TABLA, necesaria después para fichero circuits ####
write.table(comp, file = "comp.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

# Ordenar los resultados por p.valor ajustado de más a menos significativo
ranked_comp<-comp[order(comp$FDRp.value,decreasing = FALSE),]
head(ranked_comp)
#Pasar a columna los nombres de los circuitos para guardar el fichero y no perder
#esta información
name_path<-rownames(ranked_comp)
comp_ranked<-cbind(name_path,ranked_comp)
head(comp_ranked)
dim(comp_ranked)
#### ESCRIBIR ESTA TABLA DE RESULTADOS ####
write.table(comp_ranked, file = "Circuitos_dif_activ_wilcoxon_GSE164416.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)

# Computa un resumen de los resultados, nº de circuitos regulados up- y down-
pathways_summary <- get_pathways_summary(comp, metaginfo, conf = 0.05)
#View(pathways_summary)
#Pasar a columna los nombres de los circuitos para guardar el fichero y no perder
#esta información
names_sumary<-rownames(pathways_summary)
pathways_summary<-cbind(names_sumary,pathways_summary)
head(pathways_summary)
#View(pathways_summary)
#### ESCRIBIR RESUMEN RESULTADOS ####
write.table(pathways_summary, file = "Resumen_circ_desregul_GSE164416.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)

#Filtrar sólo los circuitos-efectores que pasan el umbral p.adj con FDR
comp_sig_05<-comp_ranked[comp_ranked$FDRp.value<0.05,] #Coger solo la fila con 
#adj.pvalor menor de 0.05 y todas las columnas
head(comp_sig_05)
rownames(comp_sig_05)
dim(comp_sig_05)
#157   6
comp_sig_01<-comp_ranked[comp_ranked$FDRp.value<0.01,] #Coger solo la fila con 
#adj.pvalor menor de 0.01 y todas las columnas
dim(comp_sig_01)
#29   6
#lista de circuitos-efectores significativos FDR<0.05 para hacer una aproximación
#al enriquecimiento funcional con esos efectores (script enriq_fxal_efectores_GSE164416.R)
efectors <- comp_sig_05$path_names
length(efectors)
write.table(efectors, file = "Lista_nombres_efectores_sign_wilcoxon.txt",sep = "\t", row.names = FALSE, col.names = FALSE)

###########################################################################################

##### Gráficas de resultados de activación de rutas:
#Para visualizar gráficas de resultados.
#(incluye gráfica PCA y HeatMap de los valores de activación de circuitos)

#Input data 
# Análisis PCA (orden matrix igual que orden sample_group)
ranked_path_vals<-path_vals[order(comp$FDRp.value,decreasing = FALSE),]
head(ranked_path_vals)
pca_model<-do_pca(ranked_path_vals[1:ncol(ranked_path_vals),])
pca_plot(pca_model,sample_group,legend = TRUE, main = "PCA plot", cex = 2)
#PCA multiple con varianza explicada acumulada
multiple_pca_plot(pca_model, sample_group, cex = 2, plot_variance = TRUE, comps = 3,
                  main = "Multiple PCA plot")
## Heat Map Visualización 
heatmap_plot(path_vals, group = sample_group, colors = "classic", variable_clust = TRUE, 
             labRow=rownames(path_vals), labCol=colnames(path_vals), legend = TRUE, 
             legend_xy = "topleft", main = "Heatmap")

###########################################################################################

##### Pathways comparison "rutas KEGG" pathways
#Para añadir los nodos coloreando los genes diferencialmente expresados
colors_de_hipathia<-node_color_per_de(results_physiologicals, metaginfo,
                                      sample_group,"T2D","ND", colors="hipathia")
## Visualización a través de un servidor del Report de Hipathia con todos los pathways
report_colors<-create_report(comp,metaginfo,"save_colors",node_colors = colors_de_hipathia,
                             group_by = "pathway", conf = 0.05, verbose = FALSE)
visualize_report(report_colors)

servr::daemon_stop(1)

#Para ver un pathway en concreto;
#pathway_comparison_plot(comp,metaginfo,pathway = "hsa05200",node_colors = colors_de_hipathia)

####################################################################################################

######HIPATHIA Análisis basado en funciones (UNIPROT y GO)
### Para anotación GO y Uniprot de los circuitos-efectores obtenidos de aquellos 
#desregulados en DM2:
#A partir de la variable results creada con valores de activación de rutas,
#crear dos variables nuevas para Uniprot y GO:

#### GO
GO_vals <- quantify_terms(results_physiologicals, metaginfo, "GO")
#View(assay(GO_vals))
#Utilizar GO_vals como si fuera path_vals:
### Análisis de activación diferencial de funciones con test estadístico de Wilcoxon
compGO <- do_wilcoxon(GO_vals, sample_group, g1="T2D", g2="ND", 
                    paired = FALSE, adjust = TRUE)
# Ordenar los resultados por p.valor ajustado de más a menos significativo
ranked_compGO<-compGO[order(compGO$FDRp.value,decreasing = FALSE),]
head(ranked_compGO)
#View(ranked_compGO)
write.table(ranked_compGO, file = "GO_funct_dif_wilcoxon_GSE164416.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)

# Heatmap
heatmap_plot(GO_vals, group = sample_group, colors="classic", variable_clust = TRUE)
# Análisis PCA (orden matrix igual que orden sample_group)
ranked_go_vals<-GO_vals[order(comp$FDRp.value,decreasing = FALSE),]
head(ranked_go_vals)
pca_model_GO<-do_pca(ranked_go_vals[1:ncol(ranked_go_vals),])
pca_plot(pca_model_GO,sample_group,legend = TRUE, main = "PCA plot", cex = 2)
#PCA multiple con varianza explicada acumulada
multiple_pca_plot(pca_model_GO, sample_group, cex = 2, plot_variance = TRUE, comps = 3,
                  main = "Multiple PCA plot")


#### Uniprot
Uniprot_vals <- quantify_terms(results_physiologicals, metaginfo, "uniprot")
#View(assay(Uniprot_vals))
### Análisis de activación diferencial de funciones con test estadístico de Wilcoxon
compUniprot <- do_wilcoxon(Uniprot_vals, sample_group, g1="T2D", g2="ND", 
                      paired = FALSE, adjust = TRUE)
# Ordenar los resultados por p.valor ajustado de más a menos significativo
ranked_compUniprot<-compUniprot[order(compUniprot$FDRp.value,decreasing = FALSE),]
head(ranked_compUniprot)
ranked_compUniprot$Uniprot<-rownames(ranked_compUniprot)
#View(ranked_compUniprot)
write.table(ranked_compUniprot, file = "Uniprot_dif_wilcoxon_GSE164416.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)

# Heatmap Uniprot
heatmap_plot(Uniprot_vals, group = sample_group, colors="classic", variable_clust = TRUE)
# Análisis PCA Uniprot (orden matrix igual que orden sample_group)
ranked_Uniprot_vals<-Uniprot_vals[order(compUniprot$FDRp.value,decreasing = FALSE),]
head(ranked_Uniprot_vals)
pca_model_Uni<-do_pca(ranked_Uniprot_vals[1:ncol(ranked_Uniprot_vals),])
pca_plot(pca_model_Uni,sample_group,legend = TRUE, main = "PCA plot", cex = 2)
#PCA multiple con varianza explicada acumulada
multiple_pca_plot(pca_model_Uni, sample_group, cex = 2, plot_variance = TRUE, comps = 3,
                  main = "Multiple PCA plot")





################################################################################






######################################################################
############Enriquecimiento funcional efectores circuitos#############
###########################ClusterProfilerR###########################
######################################################################

###### Leer datos de Diseño del experimento
genes <- read.table("~path/Lista_nombres_efectores_sign_wilcoxon.txt",header=FALSE,sep="\t")
head(genes)
genes<-genes$V1 #asi es un vector de cadena de caracteres con 157 circuitos (objeto tipo lista)
class(genes)
genes

#quiero sacar las proteínas efectoras para el enriquecimiento funcional
?unlist #BiocGenerics package, a partir de objeto tipo lista produce un objeto tipo vector obtenido
?strsplit #strsplit base package, split los elementos de un vector de caracteres en subcadenas de acuerdo al match
separacion <- unlist(strsplit(genes, "[:]-"))
separacion #es una lista de caracteres
class(separacion)
separacion <- data.frame(separacion) #para pasarlo a tabla
separacion
class(separacion)
#Selecciono solo las filas pares pq tienen los Symbols
seleccion <- separacion[c(FALSE,TRUE),]
seleccion
class(seleccion)
length(seleccion)
#Separo los que tienen dos o más efectores juntos (separado por espacio)
separacion2 <- unlist(strsplit(seleccion, "\\s"))
efectores<-data.frame(separacion2)
head(efectores)
dim(efectores)
# Escribir estos resultados en un fichero que se puede abrir con una hoja de calculo tipo excel, por ejemplo
write.table(efectores,file='Lista_efectores_Hipathia.txt',sep='\t', row.names=FALSE, col.names = FALSE)
#Tengo la lista de 176 efectores con los que quiero hacer enriquecimiento funcional ORA:
#Para analisis de enriquecimiento funcional GSEA necesito ID y FoldChange (geneList) 
#más adelante dice como obtener la geneList

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
#Symbol ID
#Necesito tener una lista con los nombres de los efectores:
efectores<-data.frame(genes)
dim(efectores)
head(efectores)

### GO over-representation analysis (ORA) ###
#La función enrichGO() para hacer el test de sobrerrepresentación de gene ontology
# Convert gene IDs to ENTREZID
# We will lose some genes here because not all IDs will be converted
ids<-bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
head(ids)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

#Analisis de sobrerrepresentación de GO
go_enrich <- enrichGO(gene = dedup_ids$ENTREZID, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "all",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE, pool = FALSE)
#View(as.data.frame(go_enrich))

##Barplot
library(ggplot2)
barplot(go_enrich, split="ONTOLOGY",
        drop = TRUE, showCategory = 8, 
        title = "GO Biological Pathways",
        font.size = 6, order=T)+facet_grid(ONTOLOGY~.,scale="free")
##Dotplot
dotplot(go_enrich, showCategory = 10, font.size=5, orderBy = "GeneRatio", 
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
barplot(wiki, showCategory = 15, order = T,
        title = "Enriched Pathways on WikiPathways", font.size = 8)
## Dotplot WikiPathways
dotplot(wiki, showCategory = 15, 
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
barplot(react, showCategory = 15, order = T,
        title = "Enriched Pathways on REACTOME", font.size = 6)
## Dotplot REACTOME
dotplot(react, showCategory = 20, 
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
#View(as.data.frame(kk))
#Bar Plot
barplot(kk, showCategory=14,label_format = 30, font.size = 7, 
        title = "Enriched Pathways on KEGG", order=T)
#Dot plot
#Semejante a bar plot con la capacidad de encode otros score como dot size.
dotplot(kk, showCategory=20, font.size= 7, order=T) + ggtitle("dotplot KEGG")
