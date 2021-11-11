library(gplots)
library(ggplot2)
library(RColorBrewer)

library(tidyr)

library(DEGreport)

library(gage)

library(DEGreport)
library(clusterProfiler)
library(plyr)

library(edgeR)




#############################
#DEseq2 Analysis 
#############################
library(DESeq2) #
library(qvalue)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)



#Generate countdata table from subread count table
countdata <- read.table("./Data/Subread.count", header=TRUE, row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

#Name columns
colnames(countdata) <- c("A1","A2","A3","A4","A5","A6","A7","A8",
                         "B1","B2","B3","B4","B5","B6","B7","B8")

# Convert to matrix
countdata <- as.matrix(countdata)

# Assign conditions and relevel
(genotype <- factor(c(rep("WT", 4), rep("KO", 4),rep("WT", 4), rep("KO", 4))))
genotype = relevel( genotype, "WT")
condition <- factor(c(rep("mock", 8), rep("infected", 8)))
condition = relevel( condition, "mock")

(group <- factor(c(rep("WT_mock", 4), rep("KO_mock", 4),rep("WT_inf", 4), rep("KO_inf", 4))))
group = relevel( group, "KO_mock")
group = relevel( group, "WT_inf")
group = relevel( group, "WT_mock")

#Dataframe defining the different factors for DEseq2 analysis 
(coldata <- data.frame(row.names=colnames(countdata), genotype, condition, group))

#Deseq2 dataset with a 1 factor design with 4 groups
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group)

#filter out genes where there are less than 3 samples with normalized counts greater than or equal to 3
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 3 ) >= 3
dds <- dds[idx,]

# Run the DESeq pipeline for pairwise comparisons
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

#KO vs WT in mock-infected
res_KOWT_mock <- results(dds, contrast=c("group","KO_mock","WT_mock"))    

#Infected vs mock-infected in WT and KO cells
res_Infection_in_WT <- results(dds, contrast=c("group","WT_inf","WT_mock"))
res_Infection_in_KO <- results(dds, contrast=c("group","KO_inf","KO_mock"))

# Run the DESeq pipeline with Likelihood ratio test (LRT), to idenitfy genes that are differentially regulated between any of the four conditions
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
res_LRT <- results(dds_lrt, alpha=0.01) 

#Calculate q-values
res_LRT$qvalue <- qvalue(res_LRT$pvalue)$qvalue
res_KOWT_mock$qvalue<- qvalue(res_KOWT_mock$pvalue)$qvalue
res_Infection_in_WT$qvalue<- qvalue(res_Infection_in_WT$pvalue)$qvalue
res_Infection_in_KO$qvalue<- qvalue(res_Infection_in_KO$pvalue)$qvalue

## Merge with normalized count data to generate final table
resdata <- merge(as.data.frame(res_KOWT_mock), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
rownames(resdata)<-resdata$Gene

#Add gene symbols 
resdata$symbol <-mapIds(org.Hs.eg.db,
                        keys=str_replace(resdata$Gene,pattern = ".[0-9]+$",replacement = ""), 
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resdata["SARS-CoV-2",]$symbol <- "SARS-CoV-2"

#Keep only symbol, basemean, l2fc between WT-mock and KO-mock", and q value
resdata <- resdata[,c(25,2,3,8,seq(9,24))]
names(resdata)[3] <- "l2fc_WT-mock vs KO-mock"
names(resdata)[4] <- "qvalue_WT-mock vs KO-mock"
#Add l2fc and qvalue from the WT-mock vs WT-infected and KO-mock vs KO-infected"
resdata <- cbind(resdata,as.data.frame(res_Infection_in_WT)[,c(2,7)],as.data.frame(res_Infection_in_KO)[,c(2,7)])
names(resdata)[21] <- "l2fc_WT-mock vs WT-infected"
names(resdata)[22] <- "qvalue_WT-mock vs WT-infected"
names(resdata)[23] <- "l2fc_KO-mock vs KO-infected"
names(resdata)[24] <- "qvalue_KO-mock vs KO-infected"
#reorder column
resdata <- resdata[,c(seq(1,4),seq(21,24),seq(5,20))] 
#Generate final data table
write.csv(resdata, file="diffexpr-results_forpaper.csv")

#QC of data
# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Colors for plots below
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(group))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))

png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis with rld_pca function
pdf("qc-pca.pdf")
rld_pca(rld, colors=mycols, intgroup="group", xlim=c(-7, 7), ylim=c(-7, 7))
dev.off()

# Generate Volcanoplots and Gene Ontology/KEGG pathway analysis for pairwise comparisons, using MakeGraphs function
MakeGraphs(res_KOWT_mock,dds)
MakeGraphs(res_Infection_in_WT,dds)
MakeGraphs(res_Infection_in_KO,dds)

pdf("Volcano_KOWT_mock.pdf")
Volcanoresdata(res_KOWT_mock,dds)
dev.off() 
pdf("Volcano_Infection_in_WT.pdf")
Volcanoresdata(res_Infection_in_WT,dds)
dev.off()
pdf("Volcano_Infection_in_KO.pdf")
Volcanoresdata(res_Infection_in_KO,dds)
dev.off()

#Plot individual genes
plotGene("SARS-CoV-2", dds, res_KOWT_mock)

#############################################
#Clustering of differentially expressed genes
#############################################
library(DEGreport)
library(clusterProfiler)

#select genes differentially expressed between any conditions
select_LRT <- rownames(res_LRT[ res_LRT$qvalue <= .01  & res_LRT$baseMean >15 ,],)

#Get Matrix of log normalized value, and select the differnetially expressed genes
rld_mat <- assay(rld)
Selected_rld <- rld_mat[select_LRT,]

#Generate pattern using DegPattern
clusters <- degPatterns(Selected_rld, metadata = coldata, time = "group", col='genotype')

#Add official gene symbol and ENTREZ gene id to the cluster df
clusters$df$symbol <-mapIds(org.Hs.eg.db,
                            keys=str_replace(clusters$df$genes,pattern = ".[0-9]+$",replacement = ""), 
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first") 
clusters$df$entrez <-mapIds(org.Hs.eg.db,
                            keys=str_replace(clusters$df$genes,pattern = ".[0-9]+$",replacement = ""), 
                            column="ENTREZID",
                            keytype="ENSEMBL",
                            multiVals="first") 

#Generate cluster plots (not reordered)
p <- clusters$plot +
  theme_classic() +
  theme(aspect.ratio=0.9)
p
ggsave("RNA_Cluster.pdf")

#Reorder clusters
pos.1 <- which(clusters$df$cluster ==1 )
pos.2 <- which(clusters$df$cluster ==2 )
pos.3 <- which(clusters$df$cluster ==3 )
pos.4 <- which(clusters$df$cluster ==4 )
pos.5 <- which(clusters$df$cluster ==5 )
pos.6 <- which(clusters$df$cluster ==6 )
pos.7 <- which(clusters$df$cluster ==7 )
pos.8 <- which(clusters$df$cluster ==8 )

clusters$df$cluster[pos.1] <- 2 
clusters$df$cluster[pos.2] <- 3
clusters$df$cluster[pos.3] <- 5
clusters$df$cluster[pos.4] <- 7
clusters$df$cluster[pos.5] <- 4
clusters$df$cluster[pos.6] <- 8
clusters$df$cluster[pos.7] <- 6
clusters$df$cluster[pos.8] <- 1

clusters$df <- clusters$df[order(clusters$df$cluster), ]
clusters$df["SARS.CoV.2","genes"] <- "SARS-CoV-2"

#Generate heatmap of clustered genes (reordered)

assayCluster <- assay(rld)[clusters$df$genes,]

my_palette <- colorRampPalette(c("blue",'white','red'))(n=1000)
pdf("heatmapcluster.pdf", 8, 14)
heatmap.2(assayCluster, col=my_palette,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=0.6, labRow=F, Colv=F, Rowv= F,dendrogram="none",
          main="p < .05 and log2fc >1 ")
dev.off()

#############################################
#Enrichment Analysis of clustered genes
#############################################
#Inspired from https://github.com/GoldfarbLab/H522_paper_figures

library(msigdbr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

#load Reactome and Hallmark gene sets from msig package
msig <- msigdbr(species = "Homo sapiens") %>% filter(gs_cat == "H" | (gs_cat == "C2" & gs_subcat == "CP:REACTOME"))   %>% select(gs_name,human_gene_symbol)

#Find gene sets enriched in cluster 1 to 8
enrichment.results <- tibble()
for (clusterID in c(1, 2, 3,4,5,6, 7,8)) {
    geneIds.cluster <-  clusters$df[(clusters$df$cluster ==clusterID),]
  em <- enricher(gene=geneIds.cluster$symbol,TERM2GENE=msig, minGSSize = 8, qvalueCutoff=0.05)
  #write.table(clusters$df$symbol[(clusters$df$cluster ==clusterID)],sep=" ",row.names=F,col.names=F ,file=paste("cluster_",clusterID,".txt"))
  
  if (!is.null(em)) {
    result <- em@result
    result$cluster <- clusterID
    result <- filter(result, qvalue <= 0.05)
    if (nrow(result) > 0) {
      result$genes <- apply(result, 1, function(x) {str_c(sort(unlist(str_split(x[["geneID"]],"/"))),collapse=",")} )
      enrichment.results <- rbind(enrichment.results, result)
    }
  }
}

enrichment.results$category <- ""
enrichment.results$category[str_detect(enrichment.results$Description, "^HALLMARK")] <- "Hallmark"
enrichment.results$category[str_detect(enrichment.results$Description, "^REACTOME")] <- "Reactome"
enrichment.results$category[str_detect(enrichment.results$Description, "^GO")] <- "GO_CC"
table(enrichment.results$category)

#Exemple to get enriched gene sets and genes in cluster 8
enrichment.results$Description[enrichment.results$cluster == 8]
enrichment.results$genes[enrichment.results$cluster == 8]

num.in.cluster <- data.frame(table(clusters$df$cluster))
names(num.in.cluster) <- c("cluster", "n")
enrichment.results$ratio <- enrichment.results$Count / num.in.cluster$n[enrichment.results$cluster]

# convert to matrix format
enrichment.results_wide = pivot_wider(enrichment.results, id_cols = c(Description, category), names_from = cluster, values_from=c(pvalue, p.adjust, qvalue, Count, GeneRatio, BgRatio, ratio, genes))

# The following is to make the complex heatmap of Fig 5E
#DataSet to include in the plot
include_sets <- c(
  #"REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", 
  "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
  "REACTOME_RHO_GTPASE_CYCLE",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "REACTOME_M_PHASE",
  "REACTOME_TRANSLATION",
  "REACTOME_NONSENSE_MEDIATED_DECAY_NMD", 
  "REACTOME_CELL_CYCLE_MITOTIC",  
  "REACTOME_RRNA_PROCESSING",
  "REACTOME_DNA_REPLICATION" ,
  #"REACTOME_INTERFERON_ALPHA_BETA_SIGNALING", 
  "REACTOME_NEUTROPHIL_DEGRANULATION",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_ER_QUALITY_CONTROL_COMPARTMENT_ERQC",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "REACTOME_CHOLESTEROL_BIOSYNTHESIS",
  "REACTOME_GLUCOSE_METABOLISM",
  "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_MTORC1_SIGNALING" ,
  "HALLMARK_HYPOXIA" ,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" ,
  #"HALLMARK_ESTROGEN_RESPONSE_EARLY",
  #"HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COAGULATION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" 
)

#Rename the sets
enrichment.results.pruned <- filter(enrichment.results_wide, Description %in% c(include_sets))
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_CELL_CYCLE_CHECKPOINTS")] <- "Cell cycle checkpoints"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES")] <- "Metabolism of Amino Acids"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_RHO_GTPASE_CYCLE")] <- "Signaling by Rho GTPases"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_CELL_CYCLE_CHECKPOINTS")] <- "Cell Cycle Checkpoints"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_M_PHASE")] <- "Mitosis"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_TRANSLATION")] <- "Translation"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD")] <- "Nonsense Mediated Decay"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_RRNA_PROCESSING")] <- "Ribosomal RNA processing"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_DNA_REPLICATION")] <- "DNA Replication"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_NEUTROPHIL_DEGRANULATION")] <- "Neutrophil Degranulation"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_PROTEIN_SECRETION")] <- "Protein Secretion"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")] <- "TNFa Signaling via NFkB"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC")] <- "Antigen presentation by MHC class I"    
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_ER_QUALITY_CONTROL_COMPARTMENT_ERQC")] <- "ER quality control compartment"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")] <- "Extracellular matrix organisation"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_CHOLESTEROL_BIOSYNTHESIS")] <- "Cholesterol Biosynthesis"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_GLUCOSE_METABOLISM")] <- "Glucose Metabolism"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT")] <- "TCA Cycle"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_INTERFERON_ALPHA_RESPONSE")] <- "Interferon Alpha and Beta response"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_GLYCOLYSIS")] <- "Glycolysis"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_CHOLESTEROL_HOMEOSTASIS")] <- "Cholesterol Homeostasis"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_PROTEIN_SECRETION")] <- "Protein Secretion"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_MTORC1_SIGNALING")] <- "mTORC1 Signaling"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_HYPOXIA")] <- "Hypoxia"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_INTERFERON_GAMMA_RESPONSE")] <- "Interferon Gamma response"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_FATTY_ACID_METABOLISM")] <- "Fatty Acid Metabolism"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_P53_PATHWAY")] <- "P53 Pathway"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_COAGULATION")] <- "Coagulation"
enrichment.results.pruned$Description[which(enrichment.results.pruned$Description == "HALLMARK_OXIDATIVE_PHOSPHORYLATION")] <- "Oxidative Phosphorilation"

enrichment.results.pruned <- arrange(enrichment.results.pruned, category)
em.q.values <- select(enrichment.results.pruned, grep("qvalue", colnames(enrichment.results.pruned), value=T))

em.q.values <- em.q.values[, c("qvalue_1", "qvalue_2", "qvalue_3", "qvalue_4", "qvalue_5", "qvalue_6", "qvalue_7","qvalue_8")]
rownames(em.q.values) <- as.character(enrichment.results.pruned$Description)
colnames(em.q.values) <- c(1,2,3,4,5,6,7,8)
em.q.values <- -log10(as.matrix(em.q.values))

ratio.values <- select(enrichment.results.pruned, grep("ratio", colnames(enrichment.results.pruned), value=T))
ratio.values <- ratio.values[, c("ratio_1", "ratio_2", "ratio_3", "ratio_4", "ratio_5", "ratio_6", "ratio_7","ratio_8")]

colnames(ratio.values) <- c(1,2,3,4,5,6,7,8)
ratio.values <- data.frame(ratio.values)

num.clusters <-8
max.q = 5 #max(em.q.values, na.rm=T)
min.q = min(em.q.values, na.rm=T)
#col_fun = colorRamp2(c(min.q, max.q, max.q), c("#300101", "#e41a1c", "#e41a1c"))
col_fun = colorRamp2(c(0, 2, 10), c("white",  "#e41a1c", "black"))

col.annot <- HeatmapAnnotation(cluster = anno_simple(colnames(em.q.values),
                                                     pt_gp = gpar(fontsize = 6),
                                                     height = unit(1, "mm"),
                                                     col = structure(brewer.pal(num.clusters, "Set3"), 
                                                                     names = colnames(em.q.values))),
                               show_annotation_name = F,
                               show_legend = F,
                               annotation_name_gp = gpar(fontsize = 6))


h <- Heatmap(em.q.values, 
             # labels
             column_title = "Enriched gene sets",
             column_title_gp = gpar(fontsize=7),
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 6),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize=7),
             #name = "-log10(q-value)",
             # legends
             bottom_annotation = col.annot,
             show_heatmap_legend = T,
             col = col_fun,
             heatmap_legend_param = list(color_bar = "continuous",
                                         title_gp = gpar(fontsize = 6),
                                         labels_gp = gpar(fontsize = 6),
                                         grid_width = unit(2,"mm")),
             
             cluster_rows = F,
             cluster_columns = F,
             show_row_names = T,
             show_column_names = T,
             #split = enrichment.results.pruned$category,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(fill = "white", col = "#EEEEEE"))
               grid.circle(x = x, y = y, r = sqrt(ratio.values[i, j]) * 0.1,
                           gp = gpar(fill = col_fun(em.q.values[i, j]), col = NA))
             }
)

gb_heatmap = grid.grabExpr(draw(h), height=5, width=2.5)

pdf("Enrichment.pdf", 8, 14)
h
dev.off()

##########################################################################################
#Analysis of interferon-stimulated genes and restriction factors from cluster 7-8
##########################################################################################

#Take all genes in the hallmarck gene sets "interferon alpha stimulated genes
IFNA <- as.vector(read.table("./Hallmarck_IFNA.txt", header=F))

#Genes in cluster 7 and 8 broadly link to immune response (in hallmark datasets or curated form litterature)
curated <- as.vector(c("B2M","C1S","CD47","IFITM2","IFITM3","LPAR6","PARP14","PSME1","SAMD9L",
                       "C1R","C1S","CD9","CTSB","DUSP6","FBN1","HTRA1","MMP7","THBD","TIMP1",
                       "DDX60","EIF2AK2","IFI44","IFITM1","IL15","MX1","OAS1","PLSCR1","SAMD9",
                       "APOL6","CCL2","DDX60","EIF2AK2","IFI44","IL15","MX1","NAMPT","PLSCR1","PTGS2","SAMHD1", 
                       "ATF3","ATP2B1","B4GALT5","CCL2","CCNL1","CXCL2","KLF10","NAMPT","PMEPA1","PTGS2","SLC16A6",
                       "ATP2B1","CCL2","CD55","DCBLD2","EIF2AK2","IFITM1","IL15","NAMPT","SGMS2","TPBG",
                       "IFI6","IFIT5","IFITM1","MX1","OAS1","SAMHD1",
                       "ADAR","IFIT3","IL33","TRIM22","TLR4","STAT3","CXCL5","CD38","PARP14","PTGS2","IRF1","SNCA","ALDH3A1",
                       "LCN2","MUC5AC","MMP7","TNFS15","CAST","SLFN5","NT5E","IRAK1","CXCL8","C3","TLR4"
))

#Gene list
IFN <- unique(c(IFNA[,1],curated))

#Get the same list but with Ensembl gene names
SelectIFN <- rownames(subset(resdata[select_LRT,], is.element(symbol,IFN)))

#Get average count for each condition, from the resdata count matrix
matrix.resdata <- as.matrix(resdata[,c("A1","A2","A3","A4","A5","A6","A7","A8","B1","B2","B3","B4","B5","B6","B7","B8")])
matrix.resdata <- subset(matrix.resdata, rowMeans(matrix.resdata)> 10)

resdata.avg <- as.data.frame(cbind( rowMeans(matrix.resdata[,c("A1","A2","A3","A4")]),
                                    rowMeans(matrix.resdata[,c("B1","B2","B3","B4")]),
                                    rowMeans(matrix.resdata[,c("A5","A6","A7","A8")]),
                                    rowMeans(matrix.resdata[,c("B5","B6","B7","B8")])
))
colnames(resdata.avg) <- c('WT_mock',"WT_inf",'KO_mock',"KO_inf")

#normalized to WT_mock
matrix.norm <- matrix.resdata/resdata.avg$WT_mock

#Take subset and name it with gene symbol
Normed_IFN <- matrix.norm[is.element(rownames(matrix.norm),SelectIFN),]
rownames(Normed_IFN) <- ( mapIds(org.Hs.eg.db,
                                 keys=str_replace(rownames(Normed_IFN),pattern = ".[0-9]+$",replacement = ""),
                                 column="SYMBOL",
                                 keytype="ENSEMBL",
                                 multiVals="first"))

avg.Normed_IFN_log2 <- log2(as.data.frame(cbind( rowMeans(Normed_IFN[,c("A1","A2","A3","A4")]),
                                                 rowMeans(Normed_IFN[,c("B1","B2","B3","B4")]),
                                                 rowMeans(Normed_IFN[,c("A5","A6","A7","A8")]),
                                                 rowMeans(Normed_IFN[,c("B5","B6","B7","B8")])
)))
colnames(avg.Normed_IFN_log2) <- c('WT_mock',"WT_inf",'KO_mock',"KO_inf")

my_palette2 <- colorRampPalette(c("blue4",'cornsilk','red'))(n=1000)

pdf("heatmap_IFN.pdf", 8, 14)
heatmap.2(as.matrix(avg.Normed_IFN_log2), 
          col=my_palette2, 
          scale="none", key=T, keysize=1, symkey=T, hclustfun = function(x) hclust(x, method="complete"),
          density.info="none", trace="none", Rowv = T, Colv = F, cexCol=0.5, 
          labRow=rownames(avg.Normed_IFN_log2), dendrogram="none", main=" All ",
          colsep=1:nrow(as.matrix(avg.Normed_IFN_log2)), # Add vertical grid lines
          rowsep=1:nrow(as.matrix(avg.Normed_IFN_log2)), # Add horizontal grid lines
          sepcolor = "#EEEEEE")
dev.off()

#get heatmap in an object to extract row order
dev.off()
dev.off()
hm <- heatmap.2(as.matrix(avg.Normed_IFN_log2))

#Get the qvalues between WT-mock vs KO-mock
ifn.qvalue2 <- as.data.frame(resdata[SelectIFN,c("qvalue_WT-mock vs KO-mock","symbol")])
rownames(ifn.qvalue2) <-ifn.qvalue2$symbol
colnames(ifn.qvalue2) <- c("qvalue","symbol")

ifn.qvalue <- as.matrix(ifn.qvalue2$qvalue)
rownames(ifn.qvalue) <-ifn.qvalue2$symbol
#reorder by the order in the heatmap
ifn.qvalue <- ifn.qvalue[rownames(avg.Normed_IFN_log2[rev(hm$rowInd),]),,drop=FALSE]
ifn.qvalue.log <- -log10(ifn.qvalue)

#add qvalue for legend
ifn.qvalue.log <- -log10(rbind(ifn.qvalue,0.05,0.01,0.001,0.0001,0.00001,10^-7,10^-10))

col_fun = colorRamp2(c(1, 3, 10), c("white",  "#e41a1c", "black"))

h2 <- Heatmap(ifn.qvalue.log, 
             # labels
             column_title = "Enriched gene sets",
             column_title_gp = gpar(fontsize=7),
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 6),
             column_names_rot = 0,
             row_title_gp = gpar(fontsize=7),
             #name = "-log10(q-value)",
             # legends
             show_heatmap_legend = T,
             col = col_fun,
             heatmap_legend_param = list(color_bar = "continuous",
                                         title_gp = gpar(fontsize = 6),
                                         labels_gp = gpar(fontsize = 6),
                                         grid_width = unit(2,"mm")),
             
             cluster_rows = F,
             cluster_columns = F,
             show_row_names = T,
             show_column_names = T,
             #split = enrichment.results.pruned$category,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.rect(x = x, y = y, width = width, height = height, 
                         gp = gpar(fill = "white", col = "#EEEEEE"))
               grid.circle(x = x, y = y, r =  0.01,
                           gp = gpar(fill = col_fun(ifn.qvalue.log[i,1]), col = NA))
             }
)

gb_heatmap = grid.grabExpr(draw(h), height=5, width=0.5)

pdf("p_valueforRestrictionfactorheatmap.pdf", 8, 14)
h2
dev.off()




