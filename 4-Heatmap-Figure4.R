rm(list=ls())
load('./data/geneAccNgeneSym_heatmap.RData')
load('./data/norm_tech_batch.RData')
dat_ht <- norm_tech_batch[geneAccNgeneSym_heatmap,]
hkR1 <- c('HKR1T1','HKR1T2','HKR1T3')
hkR2 <- c('HKR2T1','HKR2T2','HKR2T3')
UI <- c('UIT1','UIT2','UIT3')
LiveR1 <- c('LiveR1T1','LiveR1T2','LiveR1T3')
LiveR2 <- c('LiveR2T1','LiveR2T2','LiveR2T3')
Live3 <- c('Live3T1','Live3T2','Live3T3')
HKR1.m <- apply(dat_ht[,hkR1],1,mean)
HKR2.m <- apply(dat_ht[,hkR2],1,mean)
UI.m <- apply(dat_ht[,UI],1,mean)
LiveR1.m <- apply(dat_ht[,LiveR1],1,mean)
LiveR2.m <- apply(dat_ht[,LiveR2],1,mean)
Live3.m <- apply(dat_ht[,Live3],1,mean)

dat_mean <- data.frame(
  HKR1 = HKR1.m,
  HKR2 = HKR2.m,
  UI = UI.m,
  LiveR1 = LiveR1.m,
  LiveR2 = LiveR2.m,
  Live3 = Live3.m
)
dat_z_score <- t(apply(dat_mean,1,scale))
colnames(dat_z_score) <- colnames(dat_mean)

geneSym <- names(geneAccNgeneSym_heatmap)
names(geneSym) <- geneAccNgeneSym_heatmap
rownames(dat_z_score) <- geneSym[rownames(dat_z_score)]
source('./func-room.R')
pdf('./figures/Final-D-heatmap-zscore-meanExpr.pdf',
    width = 3.5,
    height = 10)
gplots::heatmap.2(
  x = dat_z_score,
  col = cols.gentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  lhei = layout_c$lhei,
  trace = "none",
  density.info = 'none',
  key.title ='z-score',
  key.xlab = 'z-score'
)
dev.off()
