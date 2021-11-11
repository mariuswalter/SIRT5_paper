library(gage)
library(gageData)
data(kegg.sets.hs)
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]



rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh& (abs(log2FoldChange)>1)), points(baseMean, log2FoldChange, col="green", pch=20, cex=1.5))
  #if (labelsig) {
  # require(calibrate)
  #with(subset(res, padj<thresh & (abs(log2FoldChange)>1)), textxy(baseMean, log2FoldChange, labs=symbol, cex=textcx, col=2))
}


volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1,...) {
  par(pty="s")
  with(res, plot(log2FoldChange, -log10(qvalue), pch=20, main=main, ...))
  with(subset(res, qvalue<sigthresh ), points(log2FoldChange, -log10(qvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(qvalue), pch=20, col="orange", ...))
  with(subset(res, qvalue<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(qvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
   with(subset(res, (-log10(qvalue)>3) & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(qvalue), labs=symbol, cex=textcx, ...))
  #  with(subset(res, symbol %in% gene_list ), textxy(log2FoldChange, -log10(qvalue), labs=symbol, cex=textcx, ...))
    }
  legend(legendpos, inset=c(-1,0), xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}


MakeGraphs <- function(res,dds){
  NameString <-deparse(substitute(res))
  table(res$qvalue<0.05)
  ## Order by adjusted p-value
  resOrdered <- res[order(res$qvalue), ]
  
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  rownames(resdata)<-resdata$Gene
  resdata$symbol <-mapIds(org.Hs.eg.db,
                          keys=str_replace(resdata$Gene,pattern = ".[0-9]+$",replacement = ""), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
  resdata$entrez = mapIds(org.Hs.eg.db,
                          keys=str_replace(resdata$Gene,pattern = ".[0-9]+$",replacement = ""),
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
  
  foldchanges <- resdata$log2FoldChange
  names(foldchanges) = resdata$entrez
  
  # KEGG pathway analysis
  keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=T)
  # Look at both up (greater), down (less), and statistics.
  kegname <- paste(as.character(NameString),"keggres.csv",sep="_")
  write.table( head(keggres$greater,15), sep=',', row.names=T,col.names=T ,file=kegname)
  write.table( "less",file=kegname, append = T,sep=',',row.names=T,col.names=F )
  write.table( head(keggres$less,15),file=kegname, append = T,sep=',',row.names=T,col.names=F )
  write.table( "stats",file=kegname, append = T,sep=',',row.names=T,col.names=F )              
  write.table( head(keggres$stats,15),  file=kegname, append = T, sep=',', row.names=T,col.names=F )
  
  # Gene Ontology analysis
   gobpres <- gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
   goname <- paste(as.character(NameString),"goBP.csv",sep="_")
    write.table( head(gobpres$greater,15), sep=',', row.names=T,col.names=T ,file=goname)
   write.table( "less",file=goname, append = T,sep=',',row.names=T,col.names=F )
  write.table( head(gobpres$less,15),file=goname, append = T,sep=',',row.names=T,col.names=F )
  write.table( "stats",file=goname, append = T,sep=',',row.names=T,col.names=F )              
  write.table( head(gobpres$stats,15),  file=goname, append = T, sep=',', row.names=T,col.names=F )
  
  ix <- which.min(subset(res, abs(log2FoldChange)>1)$qvalue)
  genename <- resdata[rownames(subset(res, abs(log2FoldChange)>1)[ix,]),]$symbol
  ensname <- rownames(resdata[rownames(subset(res, abs(log2FoldChange)>1)[ix,]),])
  
  png(paste(as.character(NameString),".png",sep=""), 1000, 1000, pointsize=20)
  par(mfrow=c(2,2))
  hist(res$pvalue, breaks=50, col="grey")
  maplot(subset(resdata, (!is.na(symbol)) ), main="MA Plot",thresh=0.05)
  volcanoplot(subset(resdata, (!is.na(symbol)) & (baseMean>10)), lfcthresh=1, sigthresh=0.05, textcx=.8,main=NameString,xlim=c(-4, 4),ylim=c(0, 100))
  plotCounts(dds, gene=ensname, intgroup = "group" ,main=genename)
  #plotCounts(dds, gene=ensname, intgroup = c("condition") ,main=genename)
  dev.off()              
     
}

res <- res_Infection_in_KO

plotGene <- function (GeneName,dds, res) {
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  rownames(resdata)<-resdata$Gene
  resdata$symbol <-mapIds(org.Hs.eg.db,
                          keys=str_replace(resdata$Gene,pattern = ".[0-9]+$",replacement = ""), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")  
  resdata["SARS-CoV-2",]$symbol <- "SARS-CoV-2"
  plotCounts(dds, gene=row.names(resdata[which(resdata$symbol == GeneName), ]), intgroup = c("group"),main= GeneName, transform=F)
  
  plotCounts(dds, gene=row.names(resdata[which(resdata$symbol == GeneName), ]), intgroup = c("group"),main= GeneName, returnData = T)
}

Volcanoresdata <- function(res,dds){
  NameString <-deparse(substitute(res))
  table(res$qvalue<0.05)
  ## Order by adjusted p-value
  resOrdered <- res[order(res$qvalue), ]
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  rownames(resdata)<-resdata$Gene
  resdata$symbol <-mapIds(org.Hs.eg.db,
                          keys=str_replace(resdata$Gene,pattern = ".[0-9]+$",replacement = ""), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
  resdata$entrez = mapIds(org.Hs.eg.db,
                          keys=str_replace(resdata$Gene,pattern = ".[0-9]+$",replacement = ""),
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")
  volcanoplot(subset(resdata, (!is.na(symbol)) & (baseMean>15)), lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-4, 4), ylim=c(0,100), main=NameString)
}
