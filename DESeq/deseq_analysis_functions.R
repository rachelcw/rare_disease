# function for deseq analysis


# get geneID and return symbol-gene name.
# getGeneName <- function(geneId) {
#   # Attempt to retrieve the gene symbol from the Ensembl ID
#   geneName <- AnnotationDbi::select(org.Hs.eg.db, keys = geneId, 
#                                     columns = c("SYMBOL"))
#   geneName <- geneName
#   
#   return(geneName)
#   
# }

# Function for specific gene expression

plot_gene = function(dds,gene_name,intgroup, info = ""){
  gene = plotCounts(dds, gene_name ,intgroup, returnData = T)
  p = ggplot(gene, aes(x=condition, y = count, color = condition )) + geom_point(size = 5, alpha=0.7) +
    geom_text(aes(label = row.names(gene)), size = 3, vjust = -1.3) +
    ylim(0,max(gene$count)+300) + theme_bw(base_size = 15) +
    ggtitle(paste0(gene_name,info)) 
  plot(p)
}

sig_genes = function(dds,FC = 2,pval=0.05, name = ""){
  if (name == ""){
    res <- results(dds, pAdjustMethod = "fdr" )
  }else{
    res <- results(dds, pAdjustMethod = "fdr", name  = name )
  }
  
  res = as.data.frame(res)
  res = res[!is.na(res$padj),]
  sig = res[res$padj < pval & abs(res$log2FoldChange) > log2(FC),]
  return(list(sig,res))
  
}

run_go1 = function(list_genes,universe, padj = "fdr", ont ="ALL",
                   simp = FALSE, info ="",orgdb= "org.Hs.eg.db", return_plot=F ,
                   num_cat = 20){
  
  ego2 <- enrichGO(gene         = list_genes,
                   OrgDb         =  orgdb,
                   universe =  universe,
                   keyType       = 'SYMBOL',
                   ont           = ont,
                   pAdjustMethod = padj)
  
  if (nrow(as.data.frame(ego2)) == 0) {
    print(" no res")
    return()
  }
  
  if (simp){
    ego2 = simplify(ego2)
    
  }
    
  res_go = as.data.frame(ego2)
  p = dotplot(ego2, showCategory= 10) + theme_bw(base_size = 10) + ggtitle(info)
  plot(p)
  
  if (return_plot){
    return(list(res_go,p))
  }
  
  return(res_go)
}

