# TENETPLUS_NonMultiome

## Positional Variable
$1= RNA_Counts.csv, $2= ATAC_Counts.csv, $3= UMAP.txt, $4= RNA_Pseudo.txt, $5= 16 (number of parallel jobs), $6= 1 (historyLength)

## Terminal code for running TENET_NonMul
./TENET_NonMul RNA_Counts.csv ATAC_Counts.csv UMAP.txt RNA_Pseudo.txt 16 1 human

## Example R codes for creating input files (RNA_Counts.csv, ATAC_Counts.csv, UMAP.txt, RNA_Pseudo.txt)
coembed_sub_RNA <- subset(coembed_sub,from == "RNA")

coembed_sub_ATAC <- subset(coembed_sub,from == "ATAC")

write.csv(t(as.data.frame(coembed_sub_RNA@assays$RNA@counts)),file="RNA_Counts.csv",row.names = T,quote = F)
write.csv(t(as.data.frame(coembed_sub_ATAC@assays$peaks@counts)),file="ATAC_Counts.csv",row.names = T,quote = F)
write.table(as.numeric(pseudo_rna_only),file="RNA_Pseudo.txt",row.names = F,col.names = F,quote = F)
write.table(rbind(coembed_sub_RNA@reductions$umap@cell.embeddings,coembed_sub_ATAC@reductions$umap@cell.embeddings),file="UMAP.txt",row.names = T,col.names = T,quote = F)

## Output Files for CellMatching_Script.R
CellMatching_results.csv # $1: RNA cell, $2: Matched ATAC cell, $3: UMAP distance between RNA cell and ATAC cell
CellMatching_Plots.pdf

