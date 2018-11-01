rm(list=ls())
load('./data/norm_tech_batch.RData')
pheno <- read.csv(
  './data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
  row.names = 1,
  stringsAsFactors = FALSE
)

######
# experiment 1 and experiment 2
idx.expr1 <- grep('R1', rownames(pheno))
idx.expr2 <- grep('R2', rownames(pheno))
exp12samp <- rownames(pheno)[c(idx.expr1,idx.expr2)]


#############
## HeatKilled vs Infected 
hk_raw <- exp12samp[grep('HK',exp12samp)]
hk <- pheno[hk_raw,'IDs']
li_raw <- exp12samp[grep('LI',exp12samp)]
li <- pheno[li_raw,'IDs']

hk_vs_if <- data.frame(
  cond = factor(pheno[c(hk_raw,li_raw),'Phenotype'],
                levels = c('HeatKilled','LiveInfected')
  )
)
y <- norm_tech_batch[,c(hk,li)]
pvals <-c()
fc_ <- c()
fc__ <- c()
for(i in 1:nrow(y)){
  df_ <- data.frame(x=y[i,],y=hk_vs_if$cond)
  idx.hk <- grep('HK',names(y[i,]))
  fc_[i] <- mean(y[i,-idx.hk] - y[i,idx.hk])
  fc__[i] <- mean(y[i,-idx.hk]) - mean(y[i,idx.hk])
  pvals[i] <- t.test(x~y,data=df_,paired=TRUE)$p.value
}
names(pvals) <- rownames(y)
adj.pval <- p.adjust(pvals,method = 'fdr')
hk_vs_if_paired_t_result <- data.frame(
  geneAcc = names(pvals),
  log2FC = fc_,
  FC = 2^fc_,
  p.val.paired.ttest = pvals,
  adj.pvalue.FDR = adj.pval,
  stringsAsFactors = FALSE
)

sig <- (hk_vs_if_paired_t_result$FC >= 1.3 | hk_vs_if_paired_t_result$FC <= 0.7) &
  hk_vs_if_paired_t_result$adj.pvalue.FDR <= 0.05
file.copy(
  from = '../mtb-macrophage-LD/data/gene-id-to-genesymbol-annotation.csv',
  to = './data'
)
gene.anno <- read.csv(
  './data/gene-id-to-genesymbol-annotation.csv',
  stringsAsFactors = FALSE
)
geneAcc2geneSym <- gene.anno$GeneSymbol
names(geneAcc2geneSym) <- gene.anno$GeneAcc
sig_status <- ifelse(sig,'yes','no')

geneSymbol <- geneAcc2geneSym[rownames(hk_vs_if_paired_t_result)]

out_final <- data.frame(
  hk_vs_if_paired_t_result,
  sig_status = sig_status,
  geneSymbol = geneSymbol
)
write.csv(out_final,'./data/paired-t-test-result-25Oct18.csv')
############################
############################
pdf('./figures/B-volcano.pdf',
    width = 6,
    height = 6)
plot(y=-log10(out_final$adj.pvalue.FDR),x=out_final$log2FC,
     col=ifelse(as.character(out_final$sig_status) == 'yes','red','grey50'),
     ylab = '-log10(adj-p-value)',
     xlab = 'log2(FC)'
     )
dev.off()


###########################
###########################
load('./data/norm_tech_batch.RData')
hkR1 <- c('HKR1T1','HKR1T2','HKR1T3')
hkR2 <- c('HKR2T1','HKR2T2','HKR2T3')

LiveR1 <- c('LiveR1T1','LiveR1T2','LiveR1T3')
LiveR2 <- c('LiveR2T1','LiveR2T2','LiveR2T3')

HKR1.m <- apply(norm_tech_batch[,hkR1],1,mean)
HKR2.m <- apply(norm_tech_batch[,hkR2],1,mean)
LiveR1.m <- apply(norm_tech_batch[,LiveR1],1,mean)
LiveR2.m <- apply(norm_tech_batch[,LiveR2],1,mean)

dat_z_score <- t(apply(norm_tech_batch[,c(hkR1,hkR2,LiveR1,LiveR2)],1,scale))
colnames(dat_z_score) <- c(hkR1,hkR2,LiveR1,LiveR2)
source('./func-room.R')
colsGentleman <- function() {
  library(RColorBrewer)
  hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(20)
  return(rev(hmcol))
}
mat.all <- dat_z_score 
rownames(mat.all) <- geneAcc2geneSym[rownames(mat.all)]
pdf('./figures/B-heatmap-zscore-meanExpr.pdf',
    width = 3.5,
    height = 18)
gplots::heatmap.2(
  x = mat.all,
  col = colsGentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  #lhei = layout_c$lhei,
  lhei =c(0.2, 0.05, 3, 0.25),
  breaks = seq(from=-3, to=3, length.out=21),
  trace = "none",
  density.info = 'none',
  key.title ='z-score',
  key.xlab = 'z-score'
)
dev.off()

idx.sig <- as.character(out_final$sig_status) == 'yes'
sig.genes <- rownames(out_final[idx.sig,])
sig.mat <- dat_z_score[sig.genes,]
rownames(sig.mat) <- geneAcc2geneSym[rownames(sig.mat)]
pdf('./figures/B-heatmap-zscore-meanExpr-sig-only.pdf',
    width = 3.5,
    height = 10)
gplots::heatmap.2(
  x = sig.mat,
  col = colsGentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  #lhei = layout_c$lhei,
  lhei =c(0.2, 0.05, 2.95, 0.50),
  breaks = seq(from=-3, to=3, length.out=21),
  trace = "none",
  density.info = 'none',
  key.title ='z-score',
  key.xlab = 'z-score'
)
dev.off()
idx.up <- out_final$FC >= 1.3
up.genes <- rownames(out_final[idx.sig & idx.up,])
up.mat <- dat_z_score[up.genes,]
down.genes <- rownames(out_final[idx.sig & !idx.up,])
mat <- rbind(up.mat,down.mat)
rownames(up.mat) <- geneAcc2geneSym[rownames(up.mat)]
pdf('./figures/B-heatmap-zscore-meanExpr-sig-only-UP.pdf',
    width = 3.5,
    height = 9)
gplots::heatmap.2(
  x = up.mat,
  col = colsGentleman(),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  #lhei = layout_c$lhei,
  lhei =c(0.2, 0.05, 2.75, 0.50),
  breaks = seq(from=-3, to=3, length.out=21),
  trace = "none",
  density.info = 'none',
  key.title ='z-score',
  key.xlab = 'z-score'
)
dev.off()
file.copy(from = '../mtb-macrophage-LD/data/geneAccNgeneSym_heatmap.RData',
          to = './data')
load('./data/geneAccNgeneSym_heatmap.RData')
down.genes <- rownames(out_final[idx.sig & !idx.up,])
down.mat <- dat_z_score[down.genes,]
rownames(down.mat) <- geneAcc2geneSym[rownames(down.mat)]
pdf('./figures/B-heatmap-zscore-meanExpr-sig-only-DOWN.pdf',
    width = 3.5,
    height = 8)
gplots::heatmap.2(
  x = down.mat,
  col = colsGentleman(),
  breaks = seq(from=-3, to=3, length.out=21),
  lmat = layout_c$lmat,
  lwid = layout_c$lwid,
  #lhei = layout_c$lhei,
  lhei =c(0.2, 0.05, 2.65, 0.50),
  trace = "none",
  density.info = 'none',
  key.title ='z-score',
  key.xlab = 'z-score'
)
dev.off()