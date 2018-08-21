rm(list=ls())
load('./data/norm_tech_batch.RData')
pheno <- read.csv(
  './data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
  row.names = 1,
  stringsAsFactors = FALSE
)
### annotation
gene.anno <- read.csv(
  './data/gene-id-to-genesymbol-annotation.csv',
  stringsAsFactors = FALSE
)
geneAcc2proteinAcc <- gene.anno$ProteinAcc
names(geneAcc2proteinAcc) <- gene.anno$GeneAcc
geneAcc2geneSym <- gene.anno$GeneSymbol
names(geneAcc2geneSym) <- gene.anno$GeneAcc
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
for(i in 1:nrow(y)){
  df_ <- data.frame(x=y[i,],y=hk_vs_if$cond)
  idx.hk <- grep('HK',names(y[i,]))
  fc_[i] <- mean(y[i,-idx.hk]) - mean(y[i,idx.hk])
  pvals[i] <- t.test(x~y,data=df_,paired=TRUE)$p.value
}
names(pvals) <- rownames(y)
adj.pval <- p.adjust(pvals,method = 'fdr')
hk_vs_if_paired_t_result <- data.frame(
  geneAcc = names(pvals),
  logFC = fc_,
  p.val.paired.ttest = pvals,
  adj.pvalue.FDR = adj.pval,
  stringsAsFactors = FALSE
)
sig <- abs(hk_vs_if_paired_t_result$logFC) > 0.38 & 
  hk_vs_if_paired_t_result$adj.pvalue.FDR <= 0.05

de_hk_vs_if_paired <- data.frame(
  geneSym = geneAcc2geneSym[hk_vs_if_paired_t_result$geneAcc],
  hk_vs_if_paired_t_result,
  FC = 2^hk_vs_if_paired_t_result$logFC,
  significant = ifelse(
    sig,yes = 'yes',no = 'no'
  ),
  status.Live = ifelse(
    hk_vs_if_paired_t_result$logFC > 0,
    yes = 'up',
    no = 'down'
  ),
  stringsAsFactors = FALSE
)
de_hk_vs_if_paired <- de_hk_vs_if_paired[,c(2,1,3,6,4,5,7,8)]
save(de_hk_vs_if_paired,
     file = './data/de_hk_vs_if_paired.RData')
write.csv(
  de_hk_vs_if_paired,
  './data/de_pair_ttest_HK_vs_Live_18Augut2018.csv',
  row.names = FALSE
)
####################
rm(list=ls())
load('./data/norm_tech_batch.RData')
pheno <- read.csv(
  './data/phenotype-information-mtb-lipid-droplet-TMT-31July18.csv',
  row.names = 1,
  stringsAsFactors = FALSE
)
### annotation
gene.anno <- read.csv(
  './data/gene-id-to-genesymbol-annotation.csv',
  stringsAsFactors = FALSE
)
geneAcc2proteinAcc <- gene.anno$ProteinAcc
names(geneAcc2proteinAcc) <- gene.anno$GeneAcc
geneAcc2geneSym <- gene.anno$GeneSymbol
names(geneAcc2geneSym) <- gene.anno$GeneAcc

######################
######
# experiment 3
## Uninfected vs Infected 
ui <- rownames(pheno)[grep('UI',rownames(pheno))]
ui_raw <- pheno[ui, 'IDs']
li3 <- rownames(pheno)[grep('LIR3',rownames(pheno))]
li3_raw <- pheno[li3, 'IDs']

ui_vs_li3 <- data.frame(
  cond = factor(pheno[c(ui,li3),'Phenotype'],
                levels = c('Uninfected','LiveInfected')
  )
)

y <- norm_tech_batch[,c(ui_raw,li3_raw)]
pvals <-c()
fc_ <- c()
for(i in 1:nrow(y)){
  df_ <- data.frame(x=y[i,],y=ui_vs_li3$cond)
  idx.ui <- grep('UI',names(y[i,]))
  fc_[i] <- mean(y[i,-idx.ui]) - mean(y[i,idx.ui])
  pvals[i] <- t.test(x~y,data=df_,paired=TRUE)$p.value
}
names(pvals) <- rownames(y)
adj.pval <- p.adjust(pvals,method = 'fdr')
ui_vs_if3_paired_t_result <- data.frame(
  geneAcc = names(pvals),
  logFC = fc_,
  p.val.paired.ttest = pvals,
  adj.pvalue.FDR = adj.pval,
  stringsAsFactors = FALSE
)
sig <- abs(ui_vs_if3_paired_t_result$logFC) > 0.38 & 
  ui_vs_if3_paired_t_result$adj.pvalue.FDR <= 0.05

de_ui_vs_if3_paired <- data.frame(
  geneSym = geneAcc2geneSym[ui_vs_if3_paired_t_result$geneAcc],
  ui_vs_if3_paired_t_result,
  FC = 2^ui_vs_if3_paired_t_result$logFC,
  significant = ifelse(
    sig,yes = 'yes',no = 'no'
  ),
  status.Live = ifelse(
    ui_vs_if3_paired_t_result$logFC > 0,
    yes = 'up',
    no = 'down'
  ),
  stringsAsFactors = FALSE
)
de_ui_vs_if3_paired <- de_ui_vs_if3_paired[,c(2,1,3,6,4,5,7,8)]
save(de_ui_vs_if3_paired,
     file = './data/de_ui_vs_if3_paired.RData')
write.csv(
  de_ui_vs_if3_paired,
  './data/de_pair_ttest_Uninfected_vs_Live3_18Augut2018.csv',
  row.names = FALSE
)
