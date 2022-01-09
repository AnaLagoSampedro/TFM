################################################################################
############### Construcción del fichero Circuits "Diabetes2.rds" ##################
###################### para el aprendizaje automático ##########################
################################################################################

#Generar un data.frame con 1098 filas correspondientes a todos los circuitos y 
#sus efectores en hipathia.Estos irán en una columna denominada "hipathia", 
#esta información se obtiene tras aplicar la función get_paths_data() de hipathia.
#(Se puede obtener de cualquiera de los dos ficheros GSE164416 y GTEx).

#A partir de GTEx, por ejemplo, abrir el data.frame con los 1098 subpathways que 
#hay en hipathia sacado tras aplicar la función get_paths_data(). 
#Obtenido del Script_GTEx_hipathia.R
circuits_hipathia <- readRDS("~path/circuits_hipathia.rds")
head(circuits_hipathia)
dim(circuits_hipathia)
#1098 subpathways que están en hipathia de 79 rutas fisiológicas.
colnames(circuits_hipathia) <- "hipathia" #cambiar nombre de la columna

#Ahora, generar la segunda columna del fichero circuits "in_disease":
#con el resultado de los pathways diferencialmente activados para la enfermedad,
#obtenido del objeto comp obtenido en GSE164416 tras aplicar el análisis de 
#activación diferencial de circuitos (función do_wilcoxon(); Script_HipathiaWilcoxon.R)
comp <- read.csv("~path/comp.txt", sep = "\t", header = TRUE)
head(comp)
circuits_hipathia$in_disease<-"FALSE"
circuits_hipathia$in_disease[comp$FDRp.value<=0.05]<-"TRUE"
table(circuits_hipathia$in_disease)
head(circuits_hipathia)

#Guardar este fichero circuits, como rds y feather.
write.table(circuits_hipathia, file = "Diabetes2.txt",sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(circuits_hipathia, file = "Diabetes2.rds")
#Conversión a .feather
#install.packages("feather")
library(feather)
write_feather(circuits_hipathia, path = "Diabetes2.feather")
