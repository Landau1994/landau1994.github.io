library(pheatmap)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(igraph)
library(Seurat)
library(pbmc3k.SeuratData)

mat1 <- matrix(rgamma(1000, shape = 1) * 5, ncol = 50)
#### set gene and sample name
rownames(mat1) <- paste0("Gene",1:20)
colnames(mat1) <- paste0("Sample",1:50)
#### set group
group1 <- sample(x = colnames(mat1),size = ncol(mat1)/2,replace = F)
group2 <- setdiff(colnames(mat1),group1)
#### make group1 difference from group2
mat1[,group1] <- mat1[,group1]*20

#### make another data, set group2 larger
mat2 <- matrix(rgamma(1000, shape = 1) * 5, ncol = 50)
rownames(mat2) <- paste0("Gene",21:40)
colnames(mat2) <- paste0("Sample",1:50)
mat2[,group2] <- mat2[,group2]*20
mat <- rbind(mat1,mat2)

#### add colnames data
# Generate annotations for rows and columns
annotation_col = data.frame(
  Group = factor(c(rep("group1",25),rep("group2",25))),
  CellType = factor(c(rep("CT1",25),rep("CT2",25)))
)
rownames(annotation_col) = cbind(group1,group2)
annotation_col <- annotation_col[colnames(mat),]

annotation_row = data.frame(
  GeneClass = factor(c(rep("Set1",20),rep("Set2",20)))
)
rownames(annotation_row) = paste("Gene", 1:40, sep = "")
### avoid to be overwright
expr.mat <- mat

rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
colorPal <- colorRampPalette(c("darkgreen", "yellow","red"))
tmp.col <- adjustcolor(colorPal(10), alpha.f=.8)
pheatmap(
  mat   = log10(mat+1),
  color = tmp.col,
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  border_color = NA,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  show_colnames = TRUE,
  show_rownames = TRUE,
  drop_levels   = TRUE,
  fontsize  = 8,
  main = "Pheatmap log normalized,Clustered,default color"
)


g <- random.graph.game(100,0.02)
plot.igraph(g,vertex.label="",vertex.color=tmp.col)


data("pbmc3k.final")
pbmc3k <- UpdateSeuratObject(pbmc3k.final)

length(levels(pbmc3k))
DimPlot(pbmc3k,pt.size = 3,reduction = "umap")+
  scale_color_manual(values = adjustcolor(colorPal(9), alpha.f=1))

?viridis_pal
colorPal <- colorRampPalette(c("darkgreen","#0FD64F","#F8EF42","yellow","orange","coral","red"))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
rdbu <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))
n = 265
image(
  1:n, 1, as.matrix(1:n),
  col = adjustcolor(colorPal(n), alpha.f=.7),
  xlab = "colorPal n", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
)
image(
  1:n, 1, as.matrix(1:n),
  col = adjustcolor(jet.colors(n), alpha.f=.5),
  xlab = "colorPal n", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
)
image(
  1:n, 1, as.matrix(1:n),
  col = adjustcolor(rdbu(n), alpha.f=.7),
  xlab = "colorPal n", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
)
