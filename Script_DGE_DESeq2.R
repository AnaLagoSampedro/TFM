##################################################################
################## Análisis DGE de RNA-seq #######################
#################### Proyecto GSE164416 ##########################
################# ANA MARIA LAGO SAMPEDRO ########################
##################################################################

# Establecer el directorio de trabajo (donde se guardarán los resultados)
setwd('~path/RESULTADOS')
# Para saber si estamos correctamente en el directorio de trabajo
getwd()

## ------------------------------------------------------------------------
# Instalación de los paquetes necesarios para este análisis
BiocManager::install("DESeq2")

## ------------------------------------------------------------------------
# Cargar las librerías
library(DESeq2)

## ------------------------------------------------------------------------
# Input de los datos a analizar:
## Matriz de conteos sin normalizar:  counts_matrix(T2DvsND).txt
## Matriz de diseño:  design_matrix(T2DvsND).txt 
#Es necesario crear una matriz con el diseño de la matriz de conteos;
#donde se incluyan mínimo dos columnas, la primera con el nombre de cada 
#muestra y la segunda con la condición de interés para el análisis DGE (en 
#este caso, presentar la enfermedad DM2 (T2D) o no (ND))

## ------------------------------------------------------------------------
# Definir los directorios donde se encuentran la matriz de conteos y la matriz 
# de diseño
data_folder<-"~path/FILES"

fname<-"counts_matrix(T2DvsND).txt"
fpath<-file.path(data_folder,fname)

dname<-"design_matrix(T2DvsND).txt"
dpath<-file.path(data_folder,dname)

# Leer los datos de la matriz de conteos
expreset_raw<-read.table(fpath,header=TRUE,sep="\t")
rownames(expreset_raw)<-expreset_raw$ensembl  #Asignar a cada fila el nombre del gen
expreset<-expreset_raw[,-1]  #Eliminar la primera columna, para tener una matriz 
#de conteos exclusivamente

# Leer los datos del diseño del experimento
matrix_design <- read.table(dpath,header=FALSE,sep="\t")
head(matrix_design)
names (matrix_design) = c("ID", "Condition","name")
rownames(matrix_design)<-matrix_design$ID  #Asignar a cada fila el nombre de la muestra
matrix_design<-matrix_design[,-1]  #Eliminar la primera columna (OJO que siga 
#siendo un data.frame) 

# Ya disponemos del input de datos preparado
#head(expreset)
#head(matrix_design)

################################################################################

###### ANÁLISIS DGE #####

### 1.Filtrar por genes con muy baja expresión en la mayoría de muestras:
#Con la condición siguiente; que crea una matriz con valores lógicos TRUE si la 
#condición se cumple y FALSE si no se cumple (condición valor expresión > 0.5). 
#En general, a los valores TRUE se les asignan unos y a los valores FALSE ceros,
#por lo que, al calcular la suma por filas con la función rowSums(·) lo que se 
#hace realmente es contar cuantos valores positivos se tienen por filas en la 
#matriz de valores lógicos construida anteriormente. Tras hacer esto, sólo se 
#mantienen las filas, es decir, los genes, en los que esta suma o recuento es 
#mayor o igual a 10 (al menos 10 muestras que presentan 5 o más conteos).
filtro <- expreset >=5
filtered_data <- filtro[rowSums(filtro) >= 10, ]
dim(filtered_data)
#Se trabajará con esta matriz filtrada a partir de ahora.

#Genes que permanecen para el análisis de expresión diferencial
filtrado_genes<-expreset[rownames(filtered_data),]
dim(filtrado_genes)

### 2.Expresión diferencial con DESeq2:
#Modelo probabilístico que emplea este paquete es la Binomial negativa y 
#emplea su propia estrategia de normalización. Hace Normalización, estimación 
#de la dispersión, ajuste de los datos a un modelo lineal generalizado binomial 
#negativo y comprueba la expresión diferencial mediante test paramétrico de Wald.
#Todo a través de una única función DESeq(). Lleva a cabo de manera secuencial 
#todos los paso para el análisis completo de expresión diferencial.

#Pasos;
#Crear el objeto de la clase DESeq-DataSet con los datos de conteos y la condición:
de_diabetes<-DESeqDataSetFromMatrix(countData = filtrado_genes, 
                                    colData = matrix_design, 
                                    design = ~Condition)

#Análisis de expresión diferencial
de_fit<-DESeq(de_diabetes)

#Visualizar resultados a través de la función results()
res_diabetes<-results(de_fit)

#Se puede añadir una columna con la dirección del LogFC
res_diabetes$"UP/DOWN" <- "UP"
res_diabetes$"UP/DOWN"[res_diabetes$log2FoldChange<0] <- "DOWN"

#Selección de genes Diferencialmente Expresados (DESeq2 usa por defecto p.value 0.01)
summary(res_diabetes, alpha = 0.05)

#RESULTADOS DESeq2 
#Obtener el data.frame con los resultados sin ordenar para guardar
library(dplyr) #para usar operador pipe %>%
library(tibble) #para usar función rownames_to_column()
RES_DESeq<-as.data.frame(res_diabetes) %>% rownames_to_column("nsemble")
head(RES_DESeq)
dim(RES_DESeq)
#View(RES_DESeq)
#RESULTADOS DESeq2 ordenados por logFC
RES_padj_DESeq = res_diabetes %>% as.data.frame() %>% arrange(padj)%>% rownames_to_column("nsemble")
head(RES_padj_DESeq)
dim(RES_padj_DESeq)
write.csv(x = RES_padj_DESeq, file = "DE_DESeq.csv", row.names = FALSE)

#Identificar los genes más diferencialmente expresados, se ordenan según p
significativos_0.05 <- head(res_diabetes[order(res_diabetes$padj),],n=2358)
#Lista de genes significativos mediante DESeq2
DGE_DESeq2 <- rownames(significativos_0.05)
length(DGE_DESeq2)
#2358 genes DE

#Vector de resultados significativos up-regulados:
UP_DESeq <- significativos_0.05[significativos_0.05$`UP/DOWN` == "UP",]
DGE_DESeq2_up<-rownames(UP_DESeq)
length(DGE_DESeq2_up)
#1426 genes up DE

#Vector de resultados significativos down-regulados:
DOWN_DESeq <- significativos_0.05[significativos_0.05$`UP/DOWN` == "DOWN",]
DGE_DESeq2_down<-rownames(DOWN_DESeq)
length(DGE_DESeq2_down)
#932 genes down DE

#RESULTADOS DESeq2 ordenados por logFC
RES_FC_DESeq = res_diabetes %>% as.data.frame() %>% arrange(log2FoldChange)%>% rownames_to_column("nsemble")
head(RES_FC_DESeq)
dim(RES_FC_DESeq)

### 3. Visualización gráfica de los resultados:
#Visualizar los resultados con grafica MA-Plot (plot differences versus averages for high-throughput data)
jpeg('plotMA_DESeq2.png')
plotMA(res_diabetes,alpha=0.01,main="DESeq2",ylim=c(-5.5,6.5))
dev.off()

################
# distribution of adjusted p-values (de dyplr)
jpeg("Histograma distribución p-valores ajustados.png")
hist(res_diabetes$padj, col="lightblue", main = "Adjusted p-value distribution", xlab = "adjusted p-value")
dev.off()

# distribution of non-adjusted p-values
jpeg("Histograma distribución p-valores sin ajustar.png")
hist(res_diabetes$pvalue, col="grey", main = "Non-adjusted p-value distribution", xlab = "Non-adjusted p-value")
dev.off()

################
#Volcano plot
#this plot shows the gene fold change on the x-axis against the p-value plotted on the y-axis
#to "shrink" the \(\log2\) fold changes to remove the noise associated with fold changes coming from genes with low count levels
#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 
#This helps to get more meaningful log2 fold changes for all genes independently of their expression level.
resLFC <- lfcShrink(dds = de_fit, 
                    res = res_diabetes,
                    type = "apeglm",
                    coef = 2) # corresponds to "T2D_vs_ND" comparison

head(resLFC)

# load the library if not done yet
#BiocManager::install("EnhancedVolcano")
library("EnhancedVolcano")

# The main function is named after the package
# We use the shrunken log2 fold change as noise associated with low count genes is removed 
# Name of the column in resLFC that contains the log2 fold changes
# Name of the column in resLFC that contains the p-value
jpeg("Volcano plot DGE DESeq2.png")
EnhancedVolcano(toptable = resLFC, x = "log2FoldChange", y = "padj", lab = rownames(resLFC), xlim = -3,+3, 
                title = "Volcano plot DGE DESeq2", pCutoff = 1e-05, FCcutoff = 1)
#legend=c('Not significant','Log2 fold-change (but do not pass p-value cutoff)','Pass p-value cutoff','Pass both p-value & Log2 fold change')
dev.off()

#Visualizar los resultados con grafica MA-Plot (plot differences versus averages for high-throughput data)
jpeg("plotMA_DESeq2 Shrink.png.png")
plotMA(resLFC, alpha = 0.05,main="DESeq2 (Shrink)",ylim=c(-4,5))
dev.off()
