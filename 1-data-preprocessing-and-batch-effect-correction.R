# dir.create('./data',showWarnings = FALSE)
# dir.create('./figures',showWarnings = FALSE)
# file.copy('../../mtb-macrophage-LD/data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
#           './data/')
# file.copy('../../mtb-macrophage-LD/data/final-data-abundance-no-norm-mtb-LD.csv',
#           './data/')
# file.copy('../../mtb-macrophage-LD/func-room.R',
#           './')
# file.copy('../../mtb-macrophage-LD/data/gene-id-to-genesymbol-annotation.csv',
#           './data/')
rm(list=ls())
abund_dat <- read.csv(
  './data/final-data-abundance-no-norm-mtb-LD.csv',
  row.names = 1
)

pheno <- read.csv(
  './data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
  row.names = 1
)
label_new <- as.character(pheno$IDs)
names(label_new) <- rownames(pheno)
rownames(pheno) <- label_new[rownames(pheno)]
source('./func-room.R')
raw.mat <- as.matrix(abund_dat)
colnames(raw.mat) <- label_new[colnames(raw.mat)]
raw.mat[is.na(raw.mat)] <- 1
raw.dat <- log2(raw.mat)
pdf('./figures/A-PCA-median-normalized-data-418-proteins.pdf',
    width = 5,
    height = 5)
pcaViz(
  mat = raw.dat,
  pheno = pheno,
  title = 'Median normalized data (418 proteins)'
)
dev.off()

pdf('./figures/A-HClust-median-normalized-data-418-proteins.pdf',
    width = 6,
    height = 10)
gplots::heatmap.2(
  x = raw.dat,
  col = cols.gentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  lhei = layout_c$lhei,
  trace = "none",density.info = 'none',
  key.title ='log2(abundance)',
  key.xlab = 'log2(abundance)',
  labRow = FALSE
)
dev.off()
#############
############# NA removed
raw.dat_rmNA <- as.matrix(log2(na.omit(abund_dat)))
colnames(raw.dat_rmNA) <- label_new[colnames(raw.dat_rmNA)]
pdf('./figures/A-PCA-median-normalized-data-256-proteins-NA-removed.pdf',
    width = 5,
    height = 5)
pcaViz(
  mat = raw.dat_rmNA,
  pheno = pheno,
  title = 'Median normalized data (256) | NA removed'
)
dev.off()

pdf('./figures/A-HClust-median-normalized-data-256-proteins-NA-removed.pdf',
    width = 6,
    height = 10)
gplots::heatmap.2(
  x = raw.dat_rmNA,
  col = cols.gentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  lhei = layout_c$lhei,
  trace = "none",density.info = 'none',
  key.title ='log2(abundance)',
  key.xlab = 'log2(abundance)',
  labRow = FALSE
)
dev.off()


#####################
##################### SVA | surrogate variable analysis
##################### Remove technical and batch effects
library(sva)
norm_combat <- ComBat(dat=raw.dat_rmNA, batch=pheno$DateBatch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
norm_combat2 <- ComBat(dat=norm_combat, batch=pheno$TechincalReplicate, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
norm_tech_batch <- norm_combat2
save(norm_tech_batch,
     file = './data/norm_tech_batch.RData')
pdf('./figures/A-PCA-median-normalized-data-256-proteins-NA-removed-tech-batch-effect-removed.pdf',
    width = 5,
    height = 5)
pcaViz(
  mat = norm_tech_batch,
  pheno = pheno,
  title = 'Median normalized data (256) | NA removed and\ntechnical and bacth effect removed'
)
dev.off()

pdf('./figures/A-HClust-median-normalized-data-256-proteins-NA-removed-tech-batch-effect-removed.pdf',
    width = 6,
    height = 10)
gplots::heatmap.2(
  x = norm_tech_batch,
  col = cols.gentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  lhei = layout_c$lhei,
  trace = "none",density.info = 'none',
  key.title ='log2(abundance)',
  key.xlab = 'log2(abundance)',
  labRow = FALSE
)
dev.off()
