##################################################################
############# Análisis Enriquecimiento Funcional ORA #############
#################### Proyecto GSE164416 ##########################
################# ANA MARIA LAGO SAMPEDRO ########################
##################################################################

#A continuación del Script DGE_DESeq2.R, con los resultados obtenidos:
#Este script de R describe la implementación de un análisis ORA usando 
#el paquete clusterProfiler. Más información en:
#https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# Instalar y cargar paquetes:
#BiocManager::install("clusterProfiler", version = "3.8")
#BiocManager::install("ggplot2")
#install.packages("wordcloud")
#BiocManager::install(organism, character.only = TRUE)
#BiocManager::install("ReactomePA")

library(ReactomePA)
library(clusterProfiler)
library(wordcloud)
library(ggplot2)
library(organism, character.only = TRUE)


# Anotaciones:
#Para datos de la especie Homo sapiens es necesario instalar y cargar el paquete 
#"org.Hs.eg.db". Para ver las anotaciones disponibles ir a: 
#http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
organism = "org.Hs.eg.db"

#Preparar los datos de entrada (Input):
# Leer el input a partir de los resultados DESeq2 (DE_DESeq.csv):
df = read.csv("~path/DE_DESeq.csv", header=TRUE)
# Vector con los valores de log2 fold change 
original_gene_list <- df$log2FoldChange


#Para enriquecimiento funcional términos GO:
# Nombrar el vector con los genes correspondientes
names(original_gene_list) <- df$nsemble
# Omitir cualquier valor NA 
gene_list<-na.omit(original_gene_list)
# Ordenar la lista en orden decreciente (requerido para ClusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)
# Extraer resultados significativos (padj < 0.05)UP
sig_genes_df = subset(df, padj < 0.05)
UP_genes_df = subset(sig_genes_df, UP.DOWN == "UP")
DOWN_genes_df = subset(sig_genes_df, UP.DOWN == "DOWN")
# De los resultados significativos, filtrar los log2 fold change
genesUP <- UP_genes_df$log2FoldChange
genesDOWN <- DOWN_genes_df$log2FoldChange
# Nombrar el vector con los genes correspondientes
names(genesUP) <- UP_genes_df$nsemble
names(genesDOWN) <- DOWN_genes_df$nsemble
# Omitir cualquier valor NA 
genesUP <- na.omit(genesUP)
genesDOWN <- na.omit(genesDOWN)
# filtrar on min log2fold change (log2FoldChange > 1)
genesUP <- names(genesUP)[abs(genesUP) > 1]
genesDOWN <- names(genesDOWN)[abs(genesDOWN) > 1]

#*ENRIQUECIMIENTO FUNCIONAL, ANALISIS DE SOBRERREPRESENTACION GO*
#Antes de crear el objeto enrichGO:
#Params:  
#**Ontology** Opciones: ["BP", "MF", "CC"]  
#**keyType** This is the source of the annotation (gene ids). 
#*#The options vary for each annotation. In the example of *org.Hs.eg.db*, 
#*#the options are:   
#"ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"      
#"ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "FLYBASE"      "FLYBASECG"    "FLYBASEPROT"   
#"GENENAME"     "GO"           "GOALL"        "MAP"          "ONTOLOGY"     "ONTOLOGYALL"   
#"PATH"         "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"  

## Crear el objeto enrichGO: (elegir el del grupo MF o BP o CC, se puede hacer para ALL)
#Para gene = genes (este seleccionar genesUP o genesDOWN)

###DOWN###
go_enrich <- enrichGO(gene = genesDOWN,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "all",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.1)
head(go_enrich,20)
##Barplot
barplot(go_enrich, split="ONTOLOGY",
        drop = TRUE, showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 5, order=T)+facet_grid(ONTOLOGY~.,scale="free")
##Dotplot
dotplot(go_enrich, showCategory = 9, font.size=5, orderBy = "GeneRatio", 
        split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")

###UP###
go_enrich2 <- enrichGO(gene = genesUP,
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "all",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.1)
head(go_enrich2,20)
##Barplot
barplot(go_enrich2, split="ONTOLOGY",
        drop = TRUE, showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 5, order=T)+facet_grid(ONTOLOGY~.,scale="free")
##Dotplot
dotplot(go_enrich2, showCategory = 9, font.size=5, orderBy = "GeneRatio", 
        split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")



#Para enriquecimiento funcional anotaciones KEGG:

#Para KEGG, WikiPathways, Reactome usa IDs de ENTREZ, hay que traducirlos:
# Convertir gene IDs para enrichKEGG function
head(original_gene_list)
# Se perderán algunos genes porque no todos los IDs se convertirán
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
# eliminar IDs duplicados 
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
# Crear un nuevo dataframe df2 que tenga solo los genes que son traducidos bien con bitr
df2 = df[df$nsemble %in% dedup_ids$ENSEMBL,]
# crear una nueva columna en df2 con el correspondiente ENTREZ IDs
df2$Y = dedup_ids$ENTREZID
# Crear un vector del universo de genes
kegg_gene_list <- df2$log2FoldChange
# Nombrar el vector con los ENTREZ ids
names(kegg_gene_list) <- df2$Y
# Omitir cualquier valor NA  
kegg_gene_list<-na.omit(kegg_gene_list)
# Ordenar la lista en orden decreciente (necesario para clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
# Exctraer los resultados significativos de df2
kegg_sig_genes_df = subset(df2, padj < 0.05)
kegg_UP_genes_df = subset(kegg_sig_genes_df, UP.DOWN == "UP")
kegg_DOWN_genes_df = subset(kegg_sig_genes_df, UP.DOWN == "DOWN")
# De los resultados significativos, filtrar los log2 fold change
kegg_genesUP <- kegg_UP_genes_df$log2FoldChange
kegg_genesDOWN <- kegg_DOWN_genes_df$log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genesUP) <- kegg_UP_genes_df$Y
names(kegg_genesDOWN) <- kegg_DOWN_genes_df$Y
# Omitir cualquier valor NA 
kegg_genesUP <- na.omit(kegg_genesUP)
kegg_genesDOWN <- na.omit(kegg_genesDOWN)
# filtrar on min log2fold change (log2FoldChange > 1)
kegg_genesUP <- names(kegg_genesUP)[abs(kegg_genesUP) > 1]
kegg_genesDOWN <- names(kegg_genesDOWN)[abs(kegg_genesDOWN) > 1]

##KEGG Pathway Enrichment
## Crear objeto enrichKEGG 
#**organism** KEGG Organism Code: The full list is here: https://www.genome.jp/kegg/catalog/org_list.html (need the 3 letter code). I define this as `kegg_organism` first, because it is used again below when making the pathview plots.  
#**keyType** one of 'kegg', 'ncbi-geneid', 'ncib-proteinid' or 'uniprot'.  
kegg_organism = "hsa"

###DOWN###
kk <- enrichKEGG(gene=kegg_genesDOWN, universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
dim(kk)
head(kk,38)
##Barplot
barplot(kk, 
        showCategory = 15, order = T,
        title = "Enriched Pathways KEGG",
        font.size = 8)
## Dotplot
dotplot(kk, 
        showCategory = 18, orderBy ="GeneRatio",
        title = "Enriched Pathways KEGG",
        font.size = 8)

###UP###
kk2 <- enrichKEGG(gene=kegg_genesUP, universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
dim(kk2)
head(kk2,38)
##Barplot
barplot(kk2, 
        showCategory = 20, order = T,
        title = "Enriched Pathways KEGG",
        font.size = 8)
## Dotplot
dotplot(kk2, 
        showCategory = 17, orderBy ="GeneRatio",
        title = "Enriched Pathways KEGG",
        font.size = 8)


###Análisis enriquecimiento WikiPathways
#WikiPathways es una base de datos continuamente en actualización por la comunidad de investigadores.
#Organismos que soporta; 
#get_wp_organisms()
#Para el análisis de sobrerrepresentación con este;

###DOWN###
wiki <- enrichWP(gene = kegg_genesDOWN, universe = names(kegg_gene_list), organism = "Homo sapiens")
head(wiki,20)
dim(wiki)
##Barplot WikiPathways
barplot(wiki, showCategory = 20, order = T,
        title = "Enriched Pathways on WikiPathways", font.size = 8)
## Dotplot WikiPathways
dotplot(wiki, showCategory = 20, 
        title = "Enriched Pathways on WikiPathways",
        font.size = 7, orderBy = "GeneRatio")

###UP###
wiki2 <- enrichWP(gene = kegg_genesUP, universe = names(kegg_gene_list), organism = "Homo sapiens")
head(wiki2,20)
dim(wiki2)
##Barplot WikiPathways
barplot(wiki2, showCategory = 20, order = T,
        title = "Enriched Pathways on WikiPathways", font.size = 8)
## Dotplot WikiPathways
dotplot(wiki2, showCategory = 20, 
        title = "Enriched Pathways on WikiPathways",
        font.size = 7, orderBy = "GeneRatio")
