#################################################################
############## SCRIPT NORMALIZAR CONTEOS ########################
############# OBTENER ACTIVACION PATHWAYS #######################
########################GTEx#####################################
############### ANA MARIA LAGO SAMPEDRO #########################
#################################################################
#Debido al tamaño de este fichero, este script se desarrolla en
#supercomputador.

## Instalar Bioconductor si no está instalado ya 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Instalar paquetes implementados en Bioconductor necesarios
BiocManager::install("edgeR", force = T)
BiocManager::install("limma", force = T)
BiocManager::install("hipathia", force = T)
BiocManager::install("data.table",force = T)
install.packages("feather", repos = "http://cran.us.r-project.org")

## Cargar las librerías de esos paquetes necesarios para trabajar
library(edgeR)
library(limma)
library(hipathia)
library(data.table)
library(feather)

##Leer los datos de la matriz de counts descargados de GTEx:
expreset_raw <- fread(input = "GTEx_gene_reads.gct",
               header = TRUE, 
               skip = 2,
               sep = "\t",
               data.table = FALSE,
               showProgress = TRUE)

rownames(expreset_raw)<-expreset_raw$Name  #Asignar a cada fila el nombre del gen
expreset<-expreset_raw[,-(1:2)]  #Eliminar entonces la primera y segunda columnas

################################################################################

###### Construir objeto DGE de EdgeR para la normalización de los datos
dge <- DGEList(counts = expreset)
###### Normalizar datos crudos por método TMM con calcNormFactors
tmm<-calcNormFactors(dge,method="TMM") #calcula el factor de normalización x TMM
###### Aplicar la normalización; añadir 3 counts a cada dato para después corregir 
# por el factor de normalización y a ese resultado se le hace el log2-transform
logcpm<-cpm(tmm,prior.count=3,log=TRUE)

# eliminar de los rownames el ".numero", porque Hipathia no los procesa bien
rownames(logcpm)<-gsub("\\..*", "", rownames(logcpm))

################################################################################

###### HIPATHIA ######

#Traducir los identificadores de genes de la matriz norm al específico de Hipathia
trans_data<-translate_data(logcpm,"hsa")
hhead(trans_data)
#Poner nombres de filas como primera columna 
trans_data2 <- cbind(rownames(trans_data), as.data.frame(trans_data))
hhead(trans_data2)
#Cambiar nombre de la primera columna "rownames(trans_data)" por "index" con una 
#función de data.table
setnames(trans_data2, "rownames(trans_data)", "index")
hhead(trans_data2)
#Hay que sacar el fichero .rds de la matriz de datos normalizados y logcpm transf
#de GTEx
saveRDS(trans_data2, "expreset_Hinorm_gtexV8.rds")
expreset_Hinorm_gtexV8.rds<-readRDS("expreset_Hinorm_gtexV8.rds")
#Conversión a .feather
write_feather(as.data.frame(expreset_norm_gtexV8.rds), path = "expreset_norm_gtexV8.rds.feather")



### Obtener valores de actividad de circuitos ###

# Escalar los valores a [0,1]
exp_data<-normalize_data(trans_data)
# Cargar la información de los pathways fisiológicos de humano
#Subset de pathways sin las enfermedades (fichero adjunto)
physiological_path_list<-read.table("physiological_paths.tsv",
                                    sep = "\t", stringsAsFactors = F,
                                    header = F, quote = "")

metaginfo<-load_pathways("hsa",pathways_list = physiological_path_list[,2])

# Calcular los valores de activación de rutas de todo
results_physiologicals<-hipathia(exp_data, metaginfo, decompose=FALSE,
                                 verbose=FALSE)

# Obtener la matriz de valores de activación de circuitos (path_vals) para 
#cada sujeto
path_vals_physiologicals<-get_paths_data(results_physiologicals,matrix=TRUE)
## Para poder conservar los nombres de filas en una nueva columna:
index <- rownames(path_vals_physiologicals)
path_vals_phys <- cbind(index,path_vals_physiologicals)

#Sacar el fichero .rds de la matriz de actividad de circuitos de GTEx
saveRDS(path_vals_phys, "expreset_pathways_gtexV8.rds")
expreset_pathways_gtexV8.rds<-readRDS("expreset_pathways_gtexV8.rds")
#Conversión a .feather
write_feather(as.data.frame(expreset_pathways_gtexV8.rds), path = "expreset_pathways_gtexV8.rds.feather")


################################################################################

#Para crear ficheros circuit y genes, necesaria la siguiente información de GTEx:

#guardar la info de path_vals fisiológicos que está en hipathia que me interesan, 
#son 1098 paths:
circuits_hipathia<-rownames(path_vals_physiologicals)
circuits_hipathia<-as.data.frame(circuits_hipathia)
saveRDS(circuits_hipathia, file = "circuits_hipathia.rds")

#guardar los genes que hay en GTEx
genes_GTEx<-rownames(trans_data)
genes_GTEx<-as.data.frame(genes_GTEx)
saveRDS(genes_GTEx, file="genes_GTEx.rds")

#guardar los efectores en rutas que hay en hipathia de GTEx
genes_hipathia<-as.data.frame(metaginfo$all.genes)
saveRDS(genes_hipathia, file = "genes_hipathia.rds")

