pcaViz <- function(mat,pheno,title){
  pca_norm <- prcomp(mat,scale. = TRUE)
  pca_info <- data.frame(
    PC1 = pca_norm$rotation[,1],
    PC2 = pca_norm$rotation[,2],
    sampleID = colnames(mat),
    pheno[colnames(mat),]
  )
  per.var <- round(summary(pca_norm)[[6]][2,1:2] * 100,2)
  
  library(ggplot2)
  pca_info$BiologicalReplicate <- as.character(pca_info$BiologicalReplicate)
  p <- ggplot(pca_info,
              aes(PC1,PC2)) +
    geom_point(aes(col = Phenotype, shape = TechincalReplicate),size=4,alpha=0.65) +
    geom_text(aes(label=BiologicalReplicate),size=2,color='white') +
    xlab(paste('PC1(',per.var[1],'%)',sep='')) +
    ylab(paste('PC2(',per.var[2],'%)',sep='')) +
    ggtitle(title)
  return(p)
}

cols.gentleman <- function() {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
  return(rev(hmcol))
}

layoutConfig <- function(
  vals = c(
    c(0,3),
    c(0,5),
    c(2,1),
    c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(
    0.18,
    3.0
  ),
  lhei = c(
    0.4,
    0.025,
    5,
    0.7
  )
){
  mat <- matrix(
    vals,
    nrow = nrow,
    ncol = ncol,
    byrow = TRUE
  )
  layout_c <- list(
    lmat = mat,
    lwid = lwid,
    lhei = lhei
  )
  return(layout_c)
}

layout_c <- layoutConfig(
  vals = c(c(0,3), c(0,5), c(2,1), c(0,4)),
  nrow = 4,
  ncol = 2,
  lwid = c(0.75, 3.0),
  lhei = c(0.35, 0.03, 2.75, 0.5)
)