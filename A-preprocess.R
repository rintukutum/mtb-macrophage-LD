rm(list=ls())
file.copy(from = '../mtb-macrophage-LD/data/final-data-abundance-no-norm-mtb-LD.csv',
          to = './data')
abund.dat <- read.csv(
  './data/final-data-abundance-no-norm-mtb-LD.csv',
  row.names = 1,
  stringsAsFactors = FALSE
  )
file.copy(from = '../mtb-macrophage-LD/data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
          to = './data')

pheno <- read.csv(
  './data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
  row.names = 1,
  stringsAsFactors = FALSE
)
idx.status <- pheno$ValidationSet == 'yes'
raw.dat2 <- abund.dat[,!idx.status]
pheno.2 <- pheno[!idx.status,]
label_new <- as.character(pheno.2$IDs)
names(label_new) <- rownames(pheno.2)
rownames(pheno.2) <- label_new[rownames(pheno.2)]
file.copy(
  from = '../mtb-macrophage-LD/func-room.R',
  to = './')
source('./func-room.R')
raw.mat <- as.matrix(na.omit(raw.dat2))
colnames(raw.mat) <- label_new[colnames(raw.mat)]
raw.mat[is.na(raw.mat)] <- 1
raw.dat <- log2(raw.mat)
pdf('./figures/A-HClust-median-normalized-data-258-proteins.pdf',
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

library(sva)
norm_combat <- ComBat(dat=raw.dat, batch=pheno.2$DateBatch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
norm_combat2 <- ComBat(dat=norm_combat, batch=pheno.2$TechincalReplicate, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
norm_tech_batch <- norm_combat2
dir.create('./data',showWarnings = FALSE)
save(norm_tech_batch,
     file = './data/norm_tech_batch.RData')


pdf('./figures/A-HClust-median-normalized-data-258-proteins-NA-removed-tech-batch-effect-removed.pdf',
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

