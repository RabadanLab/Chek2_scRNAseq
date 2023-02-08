library(GSEABase)
library(AUCell)
library(effsize)

mat <- as.matrix(t(read.table(
  sep=',',
  header=TRUE,
  check.names=FALSE,
  "/Volumes/zjf/CHEK2/GSM3828672.X.csv",
  row.names=1
)))

#MSigDB <- getGmt("h_c2_c6.all.v7.4.symbols.gmt")
MSigDB <- getGmt("Geneset.gmt")

Genesets <- matrix("",nrow=length(MSigDB),ncol=3)
colnames(Genesets) <- c("setName","Source","geneIds")
for(i in 1:length(MSigDB)){
  Genesets[i,"setName"] <- setName(MSigDB[[i]])
  Genesets[i,"geneIds"] <- paste0(geneIds(MSigDB[[i]]),collapse = ", ")
}
write.csv(Genesets,"Genesets_for_Supp.csv",row.names = F)

cells_rankings <- AUCell_buildRankings(mat, nCores=6, plotStats=FALSE)

cells_AUC <- AUCell_calcAUC(MSigDB, cells_rankings,nCores = 6)
#save(cells_AUC, file=paste0("results/v02/B/interim_data/",patient,".cells_AUC.RData"))
write.csv(getAUC(cells_AUC), file="/Volumes/zjf/CHEK2/GSM3828672.cells_AUC.Geneset.csv")

##################################################################################
for(LABEL in c("Macrophage", "Malignant","Oligodendrocyte")){
  NES_mat <- data.matrix(read.csv("/Volumes/zjf/CHEK2/GSM3828672.cells_AUC.csv",
                                  header = TRUE,
                                  row.names = 1,
                                  check.names = FALSE))
  metadata <- read.table("/Users/junfeizhao/Downloads/GSE131928_RAW/IDHwtGBM.tSNE.SS2.txt",
                         sep="\t",header = T)[-1,]
  metadata <- metadata[metadata$LABELS==LABEL,]
  
  NES_mat <- NES_mat[,metadata$NAME]
  
  CHEK2_cells <- unlist(read.csv("CHEK2_cells.csv",check.names = F,header = F))
  
  CHEK2_cells <- intersect(CHEK2_cells,colnames(NES_mat))
  Non_CHEK2_cells <- setdiff(colnames(NES_mat),CHEK2_cells)
  
  All_pathways <- rownames(NES_mat)
  result <- matrix(NA,nrow=nrow(NES_mat),ncol=5,
                   dimnames=list(All_pathways,c("mean_CHEK2","mean_Non_CHEK2","effsize","pvalue","adj.p")))
  for(pathway in All_pathways){
    tmp1 <- as.numeric(NES_mat[pathway,CHEK2_cells])
    tmp2 <- as.numeric(NES_mat[pathway,Non_CHEK2_cells])
    result[pathway,"mean_CHEK2"] <- mean(tmp1,na.rm=TRUE)
    result[pathway,"mean_Non_CHEK2"] <- mean(tmp2,na.rm=TRUE)
    result[pathway,"effsize"] <- cohen.d(tmp1,tmp2)$estimate
    result[pathway,"pvalue"] <- wilcox.test(tmp1,tmp2)$p.value
  }
  result[,"adj.p"] <- p.adjust(result[,"pvalue"],method="fdr")
  write.csv(result,paste0("AUCell_DE_analysis.",LABEL,".csv"))
}

##################################################################

NES_mat <- data.matrix(read.csv("/Volumes/zjf/CHEK2/GSM3828672.cells_AUC.csv",
                                header = TRUE,
                                row.names = 1,
                                check.names = FALSE))
metadata <- read.table("/Volumes/zjf/CHEK2/GSE131928_RAW/IDHwtGBM.tSNE.SS2.txt",
                       sep="\t",header = T)[-1,]
metadata <- metadata[metadata$LABELS=="Macrophage",]

NES_mat <- NES_mat[,metadata$NAME]

CHEK2_cells <- unlist(read.csv("CHEK2_cells.csv",check.names = F,header = F))

CHEK2_cells <- intersect(CHEK2_cells,colnames(NES_mat))
Non_CHEK2_cells <- setdiff(colnames(NES_mat),CHEK2_cells)

NES_mat <- NES_mat[,c(CHEK2_cells,Non_CHEK2_cells)]

res <- read.csv("AUCell_DE_analysis.Macrophage.csv",row.names = 1)
res <- res[order(res$adj.p,decreasing = F),]
res <- res[which(res$adj.p < 0.05 & res$effsize > 0.5),]

library(pheatmap)

pdf("AUCell.Macrophage.pdf",height = 5,width = 8)

pheatmap(NES_mat[rownames(res),],
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = length(CHEK2_cells),
         scale = "row")

dev.off()



