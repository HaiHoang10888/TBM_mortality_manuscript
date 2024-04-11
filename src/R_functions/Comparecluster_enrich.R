## Deseq2
res_p <- res_DE
upgene <- subset(res_p,log2FoldChange>0)
downgene <- subset(res_p,log2FoldChange<0)
nrow(upgene)
nrow(downgene)

#######################################
up_gene = rownames(upgene)
down_gene =rownames(downgene)
Up_gene_ID <- bitr(up_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Down_gene_ID <- bitr(down_gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


genelist <- list("Upregulated" = up_gene,
                 "Downregulated" = down_gene)
ck <- compareCluster(geneCluster = genelist, keyType = "ENSEMBL",fun = "enrichGO", OrgDb="org.Hs.eg.db",ont = "BP",readable = T)
ck <- ck_old <- simplify(ck)
# jpeg("bootstrap/Figure/Comparecluste_136_DE_gene.jpeg",units = "in", width = 7, height = 5, res = 300)
# dotplot(ck,showCategory=15,font.size = 10,group=F)
# dev.off()
save(ck, file = "bootstrap/Comparecluste_136_DE_gene.RData")
cnetplot(ck)


## ============================ old 125 DE gene==============================
res_DE_dis_125 <- res_DE_dis[Common_125,]
res_p <- res_DE_dis_125
upgene <- subset(res_p,log2FoldChange>0)
downgene <- subset(res_p,log2FoldChange<0)
nrow(upgene)
nrow(downgene)

#######################################
up_gene = rownames(upgene)
down_gene =rownames(downgene)
Up_gene_ID <- bitr(up_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
Down_gene_ID <- bitr(down_gene, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


genelist <- list("Upregulated" = up_gene,
                 "Downregulated" = down_gene)
ck <- compareCluster(geneCluster = genelist, keyType = 'ENSEMBL',fun = "enrichGO", OrgDb="org.Hs.eg.db",ont = "BP",readable = T)
ck <- ck_new <- simplify(ck)
jpeg("bootstrap/Figure/Comparecluste_125_DE_gene_old_analysis.jpeg",units = "in", width = 7, height = 5, res = 300)
dotplot(ck,showCategory=15,font.size = 10,group=F)
dev.off()
save(ck, file = "bootstrap/Comparecluste_125_DE_gene_old_analysis.RData")

##=================Venn diagram old vs new
all_genes <- list("New analysis" = ck_new@compareClusterResult$ID,
                  "Old analysis" = ck_old@compareClusterResult$ID)
jpeg("bootstrap/Figure/Venn_New_vs_Old_pathway_approach1.jpeg",units = "in", width = 5, height = 5, res = 300)
create_venn_data(all_genes,main="Repeated pathways in New and Old analysis",category.names = names(all_genes))
dev.off()







ck <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
ck <-setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

####################
