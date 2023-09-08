rm(list=ls())
cat("\f")

pdf("./CellMatching_Plots.pdf",width = 8,height = 5)
#### *library*
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(progress)
library(cowplot)
library(tictoc)
library(gridExtra)

#CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
#setwd(CURRENT_WORKING_DIR)
getwd()

#################### 1. Make Coembed pseudotime
##### 1) Load Data
RNA_Counts <- read.csv("./RNA_Counts.csv",row.names = 1,header = TRUE,check.names = FALSE)
ATAC_Counts <- read.csv("./ATAC_Counts.csv",row.names = 1,header = TRUE,check.names = FALSE)
UMAP <- read.table("./UMAP.txt")
UMAP$orig <- ifelse(rownames(UMAP) %in% rownames(RNA_Counts),"RNA","ATAC")

# Plot1
ggplot()+
  geom_point(data=UMAP,aes(x=UMAP_1,y=UMAP_2,col=orig))+
  ggtitle("UMAP coordinate")+
  labs(subtitle = paste0("Total RNA cell number: ",sum(UMAP$orig=="RNA"),"\n",
                         "Total ATAC cell number: ",sum(UMAP$orig=="ATAC")))+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))

UMAP_RNA <- subset(UMAP,orig=="RNA")[,c("UMAP_1","UMAP_2")]
UMAP_ATAC <- subset(UMAP,orig=="ATAC")[,c("UMAP_1","UMAP_2")]
dim(UMAP_RNA) ; dim(UMAP_ATAC)

UMAP_RNA <- UMAP_RNA[rownames(RNA_Counts),]
UMAP_ATAC <- UMAP_ATAC[rownames(ATAC_Counts),]
table(rownames(UMAP_RNA) == rownames(RNA_Counts)) ; table(rownames(UMAP_ATAC) == rownames(ATAC_Counts))

#################### 2. Mapping with UMAP distance
##### 1) Calc distance
dist_rna <- data.frame(matrix(nrow=nrow(RNA_Counts),ncol=2))
colnames(dist_rna) <- c("ATAC_cell_id","dist")
rownames(dist_rna) <- rownames(RNA_Counts)

tic("Run time")
pb <- progress_bar$new(total = nrow(RNA_Counts))
for (j in 1:nrow(RNA_Counts)){
  RNA_umap_tmp <- UMAP_RNA[j,]
  dist_tmp <- c()
  for (k in 1:nrow(ATAC_Counts)){
    ATAC_umap_tmp <- UMAP_ATAC[k,]
    dist_tmp <- c(dist_tmp,dist(rbind(RNA_umap_tmp,ATAC_umap_tmp), method = "euclidian"))
  }
  idx=which.min(dist_tmp)
  dist_rna[j,] <- c(rownames(ATAC_Counts)[idx],dist_tmp[idx])
  pb$tick()
}
toc()
dist_rna$dist <- as.numeric(dist_rna$dist)
rm(RNA_umap_tmp,ATAC_umap_tmp,dist_tmp,idx,j,k,pb)
write.csv(dist_rna,paste0("./CellMatching_results.csv"),row.names = T,quote = F)

##### 2) Plot
RNA_Pseudo <- scan("./RNA_Pseudo.txt")
order_RNA <- order(RNA_Pseudo)
head(RNA_Pseudo[order_RNA])
head(dist_rna[order_RNA,])

# Plot2
ggplot()+
  geom_point(data=UMAP_RNA,aes(x=UMAP_1,UMAP_2,col=RNA_Pseudo))+
  scale_color_gradient(low="blue",high="yellow")+
  ggtitle("Pseudotime of RNA cells")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))

# Plot3
ggplot()+
  geom_point(data=dist_rna[order_RNA,],aes(x=1:nrow(dist_rna),y=dist),alpha=0.5)+
  xlab("pseudotime")+
  ggtitle("Distance of each cell")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))

# Plot4
cell_rna <- unique(rownames(dist_rna))
cell_atac <- unique(dist_rna$ATAC_cell_id)

grid.arrange(ncol=2,
ggplot()+
  geom_point(data=UMAP_RNA,aes(x=UMAP_1,y=UMAP_2,colour="Unselected"))+
  geom_point(data=UMAP_RNA[which(rownames(UMAP_RNA) %in% cell_rna),],
             aes(x=UMAP_1,y=UMAP_2,colour="RNA"))+
  scale_colour_manual("", 
                      breaks = c("Unselected", "RNA"),
                      values = c("grey", "light coral")) +
  ggtitle("RNA")+
  labs(subtitle = paste0("Selected RNA cell number: ",length(cell_rna)," / ",nrow(RNA_Counts)))+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))
,
ggplot()+
  geom_point(data=UMAP_ATAC,aes(x=UMAP_1,y=UMAP_2,colour="Unselected"))+
  geom_point(data=UMAP_ATAC[which(rownames(UMAP_ATAC) %in% cell_atac),],
             aes(x=UMAP_1,y=UMAP_2,colour="ATAC"))+
  scale_colour_manual("", 
                      breaks = c("Unselected", "ATAC"),
                      values = c("grey", "steel blue")) +
  ggtitle("ATAC")+
  labs(subtitle = paste0("Selected ATAC cell number: ",length(cell_atac)," / ",nrow(ATAC_Counts)))+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))
)

# Plot5
ggplot()+
  geom_point(data=UMAP,aes(x=UMAP_1,y=UMAP_2,colour="Unselected"))+
  geom_point(data=UMAP_RNA[which(rownames(UMAP_RNA) %in% cell_rna),],
             aes(x=UMAP_1,y=UMAP_2,colour="RNA"))+
  geom_point(data=UMAP_ATAC[which(rownames(UMAP_ATAC) %in% cell_atac),],
             aes(x=UMAP_1,y=UMAP_2,colour="ATAC"))+
  scale_colour_manual("", 
                      breaks = c("Unselected", "RNA", "ATAC"),
                      values = c("grey", "light coral", "steel blue")) +
  ggtitle("RNA & ATAC")+
  theme(plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))
dev.off()


#################### 4. Make TENET input files
## 1) Merge
# RNA
MatrixRNA <- RNA_Counts[rownames(dist_rna),]
dim(MatrixRNA) ; MatrixRNA[1:3,1:3]

# ATAC
MatrixATAC <- ATAC_Counts[dist_rna$ATAC_cell_id,]
dim(MatrixATAC) ; MatrixATAC[1:3,1:3]

merged_matrix <- cbind(MatrixRNA,MatrixATAC)
dim(merged_matrix) ; merged_matrix[1:3,1:3]

table(rowSums(merged_matrix)==0)
table(colSums(merged_matrix)==0)
merged_matrix <- merged_matrix[,colSums(merged_matrix)!=0]
table(colSums(merged_matrix)==0) ; dim(merged_matrix)

### 2) Save Tenet input file
write.csv(merged_matrix, file="MergedMatrix.csv",row.names = T,quote = F)
write.table(rep(1,nrow(merged_matrix)),file="CellSelect.txt",row.names = F,col.names = F,quote = F)

rm(list=ls())
