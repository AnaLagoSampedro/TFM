################################################################################
############### Construcción del fichero genes "genesANA.rds" ##################
###################### para el aprendizaje automático ##########################
################################################################################

# Establecer el directorio de trabajo (donde se guardarán los resultados)
setwd('~path/RESULTADOS')

## ------------------------------------------------------------------------

### Construcción del fichero genes ###
#Necesario generar otro data.frame que contenga tres columnas, la primera con 
#todos los identificadores de entrez, la segunda contiene solo aquellos genes 
#que aparecen en hipathia y la tercera los que son genes diana de fármacos
#aprobados en DrugBank filtrados (approved_targets).

#Cargar librerías necesarias:
pacman::p_load("biomaRt", "org.Hs.eg.db", "magrittr", 
               "AnnotationDbi", "regexPipes","here")
#Antes; 
#Filtrar genes diana de fármacos aprobados de DrugBank;
#Cargar base de datos descargada de DrugBank para filtrar aquellos genes diana
#aprobados con acción farmacológica conocida
data_folder <- "~path/FILES"
fname <- "drugbank_drug-bindings_v5.1.8.tsv"
fpath <-file.path(data_folder,fname)
drugbank_alltar <- read.delim(file = fpath, sep = "\t")
#Seleccionar aquellos fármacos aprobados y que no contengan el tag "Withdrawn"
drugbank_approved <- drugbank_alltar[-(regexPipes::grep(drugbank_alltar$groups,
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
approved_targets<-unique(drugbank_app_action$entrez_id)
#View(as.data.frame(approved_targets))


#Para la primera columna del fichero "genesANA.rds", necesarios todos los ID de 
#genes de entrez. Esto se había obtenido de GTEx tras la función translate_data() 
#y se obtenían 23664 genes (fichero obtenido en el script Script_GTEx_hipathia.R)
genes <- readRDS("~path/genes_GTEx.rds")
head(genes)
dim(genes)
#View(genes)
#Generar el data.frame con la información necesaria para el fichero genes.rds 
#para la primera columna entrezs:
entrezs<-genes$genes_GTEx
genesANA<-as.data.frame(entrezs)
#View(genesANA)
#para la segunda columna "in_hipathia":
#con la info de los genes de rutas fisiológicas que aparecen en hipathia en la 
#opción metaginfo$all.genes, obtenido del proyecto GTEx, sale todo lo que hay en 
#hipathia (fichero obtenido en el script Script_GTEx_hipathia.R)
in_hipathia <- readRDS("~path/genes_hipathia.rds")
head(in_hipathia)
dim(in_hipathia)
genesANA$in_hipathia <- as.integer(genesANA$entrezs %in% in_hipathia$`metaginfo$all.genes`)
genesANA$in_hipathia[genesANA$in_hipathia==1]<-"TRUE"
genesANA$in_hipathia[genesANA$in_hipathia==0]<-"FALSE"
#para la tercera columna "approved_targets":
#con la info de los genes diana de fármacos aprobados en la lista approved_targets 
head(approved_targets)
dim(as.data.frame(approved_targets))
genesANA$approved_targets <- as.integer(genesANA$entrezs %in% approved_targets)
genesANA$approved_targets[genesANA$approved_targets==1]<-"TRUE"
genesANA$approved_targets[genesANA$approved_targets==0]<-"FALSE"
#Visualizar los resultados
table(genesANA$in_hipathia)
table(genesANA$approved_targets)
#View(genesANA)
#Guardar este fichero genes, como rds y feather.
write.table(genesANA, file = "genesANA.txt",sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(genesANA, file = "genesANA.rds")
#Conversión a .feather
#install.packages("feather", repos = "http://cran.us.r-project.org")
library(feather)
write_feather(genesANA, path = "genesANA.feather")

################################################################################