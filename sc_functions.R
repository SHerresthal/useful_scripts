GeneTypeExpr <- function(input=expr){
  genes <- rownames(input)
  IDX <- match(genes, genetypes$SYMBOL)
  genetypeexpr <- data.frame(SYMBOL=genes,
                             TYPE=genetypes$TYPE[IDX],
                             MEAN=rowMeans(input),
                             SUM=rowSums(input))
  ggplot(genetypeexpr,aes(x=TYPE,y=SUM))+
    geom_jitter(height = 0, width = 0.1)+
    scale_y_log10()+
    ylab("Sum over all cells") +
    xlab("")+
    theme_classic() +
    coord_flip()+
    theme(text = element_text(size=12),legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))
}


highestGenes <- function(input=expr,
                         numGenes=50){
  tmp <- input
  tmp <- tmp[order(rowSums(tmp), decreasing = T),]
  tmp <- tmp[1:numGenes,]
  tmp <- melt(t(tmp))
  colnames(tmp)<- c("sample","gene","value")
  tmp$gene <- factor(tmp$gene, levels = rev(unique(tmp$gene)))
  ggplot(tmp, aes(x = tmp$gene, y = value)) +
    geom_boxplot()+
    scale_y_continuous()+
    xlab("Gene")+
    ylab("Raw UMI Counts")+
    ggtitle(paste("Counts of", numGenes, "highest expressed genes")) +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(size=8, angle = 90, hjust = 1),
          plot.title = element_text(size = 8, face = "bold"))
}

GSEA <- function(object,
                 condition_up,
                 condition_down,
                 cluster,
                 clustering,
                 condition ="status",
                 top=50,
                 GeneSets =c("GO","KEGG","DO","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                 GOntology = "BP",
                 pCorrection = "bonferroni", # choose the p-value adjustment method
                 pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                 qvalueCutoff = 0.05 # set the q-value cutoff (FDR corrected)
){
  Idents(object = object) <- clustering
  tmp <- SubsetData(object = object,ident.use = cluster)
  present_genes <- as.matrix(GetAssayData(object = tmp, slot = 'counts'))
  present_genes <- rownames(present_genes[apply(present_genes, 1, function (x) {sum(x >= 1)})>3,])
  present_genes_entrez <- bitr(present_genes, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  Idents(object = tmp) <- condition
  markers <- FindAllMarkers(object = tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers %>% group_by(cluster) %>% top_n(n = top, wt = avg_logFC) -> top
  top_up <- top[top$cluster==condition_up,]$gene
  entrez_up <- bitr(top_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  top_down <- top[top$cluster==condition_down,]$gene
  entrez_down <- bitr(top_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  OrgDb = org.Hs.eg.db
  results <- list()
  # GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                           universe = present_genes_entrez,
                                           OrgDb = OrgDb,
                                           ont = GOntology,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff  = qvalueCutoff,
                                           readable      = T))
    if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
                                             universe = present_genes_entrez,
                                             OrgDb = OrgDb,
                                             ont = GOntology,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff  = qvalueCutoff,
                                             readable      = T))
    if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO enrichment for genes downregulated in cluster ",cluster,sep="")}
  }
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    org = "hsa"
    results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up,
                                               organism = org,
                                               universe = present_genes_entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down,
                                                 organism = org,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG enrichment for genes downregulated in cluster ",cluster,sep="")}
  }
  # DO enrichment
  if("DO" %in% GeneSets){
    print("Performing Disease Ontology enrichment")
    results$DOup <- as.data.frame(enrichDO(gene = entrez_up,
                                           universe = present_genes_entrez,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff = qvalueCutoff,
                                           minGSSize     = 5,
                                           maxGSSize     = 500,
                                           readable=TRUE))
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$DOdown <- as.data.frame(enrichDO(gene = entrez_down,
                                             universe = present_genes_entrez,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff = qvalueCutoff,
                                             minGSSize     = 5,
                                             maxGSSize     = 500,
                                             readable=TRUE))
    if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in cluster ",cluster,sep="")}
  }
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    results$HALLMARKup <- as.data.frame(enricher(entrez_up,
                                                 TERM2GENE=hallmark_genes,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$HALLMARKdown <- as.data.frame(enricher(entrez_down,
                                                   TERM2GENE=hallmark_genes,
                                                   universe = present_genes_entrez,
                                                   pAdjustMethod = pCorrection,
                                                   pvalueCutoff  = pvalueCutoff,
                                                   qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK enrichment for genes downregulated in cluster ",cluster,sep="")}
  }
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    print("Performing Cannonical Pathway (C2) enrichment")
    results$cannonicalPathwaysup <- as.data.frame(enricher(entrez_up,
                                                           TERM2GENE=cannonicalPathway_genes,
                                                           universe = present_genes_entrez,
                                                           pAdjustMethod = pCorrection,
                                                           pvalueCutoff  = pvalueCutoff,
                                                           qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysup)>0){results$cannonicalPathwaysup$Enrichment <- paste("Cannonical pathway enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$cannonicalPathwaysdown <- as.data.frame(enricher(entrez_down,
                                                             TERM2GENE=cannonicalPathway_genes,
                                                             universe = present_genes_entrez,
                                                             pAdjustMethod = pCorrection,
                                                             pvalueCutoff  = pvalueCutoff,
                                                             qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysdown)>0){results$cannonicalPathwaysdown$Enrichment <- paste("Cannonical pathway enrichment for genes downregulated in cluster ",cluster,sep="")}
  }
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    results$Motifup <- as.data.frame(enricher(entrez_up,
                                              TERM2GENE=motifs,
                                              universe = present_genes_entrez,
                                              pAdjustMethod = pCorrection,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$Motifdown <- as.data.frame(enricher(entrez_down,
                                                TERM2GENE=motifs,
                                                universe = present_genes_entrez,
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in cluster",cluster,sep="")}
  }
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing immunesignature enrichment")
    results$ImmSigup <- as.data.frame(enricher(entrez_up,
                                               TERM2GENE=immuno_genes,
                                               universe = present_genes_entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in cluster ",cluster,sep="")}
    results$ImmSigdown <- as.data.frame(enricher(entrez_down,
                                                 TERM2GENE=immuno_genes,
                                                 universe = present_genes_entrez,
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in cluster ",cluster,sep="")}
  }
  return(results)
}

dotplotGSEA <- function(x,
                        show=25,
                        font.size=10,
                        title.size=10,
                        title.width=100,
                        order="count"){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    x <- if(nrow(x)>show){x[c(1:show),]}else{x}
    if(order=="padj"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- ifelse(nchar(x$Description)<100,
                              paste(substr(x$Description, 1, 100),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = unique(x$Description))
    }
    if(order=="count"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$Description <- ifelse(nchar(x$Description)>60,
                              paste(substr(x$Description, 1, 60),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = unique(x$Description))
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
    }
    ggplot(x, aes(x = GeneRatio, y = Description, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      scale_colour_gradientn(colours=c('red',
                                       'orange',
                                       'darkblue',
                                       'darkblue'),
                             limits=c(0,1),
                             values   = c(0,0.05,0.2,0.5,1),
                             breaks   = c(0.05,0.2,1),
                             labels = format(c(0.05,0.2,1))) +
      ylab(NULL) +
      ggtitle(paste(strwrap(unique(x$Enrichment), width=title.width), collapse = "\n"))+
      theme_bw() +
      theme(text = element_text(size=font.size),
            plot.title = element_text(size=title.size))
  }
}