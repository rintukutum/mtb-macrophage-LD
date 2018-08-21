rm(list=ls())
load('./data/de_hk_vs_if_paired.RData')

sig_hk<- de_hk_vs_if_paired$significant == 'yes'

idx.up_hk_vs_if <- which(de_hk_vs_if_paired[sig_hk,'status.Live'] == 'up')
idx.down_hk_vs_if <- which(de_hk_vs_if_paired[sig_hk,'status.Live'] == 'down')

genSym.up.hk_vs_if <- de_hk_vs_if_paired$geneSym[sig_hk][idx.up_hk_vs_if]
genSym.down.hk_vs_if <- de_hk_vs_if_paired$geneSym[sig_hk][idx.down_hk_vs_if]

load('./data/de_ui_vs_if3_paired.RData')

sig_ui<- de_ui_vs_if3_paired$significant == 'yes'

idx.up_ui_vs_if <- which(de_ui_vs_if3_paired[sig_ui,'status.Live'] == 'up')
idx.down_ui_vs_if <- which(de_ui_vs_if3_paired[sig_ui,'status.Live'] == 'down')

genSym.up.ui_vs_if <- de_ui_vs_if3_paired$geneSym[sig_ui][idx.up_ui_vs_if]
genSym.down.ui_vs_if <- de_ui_vs_if3_paired$geneSym[sig_ui][idx.down_ui_vs_if]

elements_down <- c(
  genSym.down.hk_vs_if,
  genSym.down.ui_vs_if
)
# HK down
length(genSym.down.hk_vs_if)
# UI down
length(genSym.down.ui_vs_if)
# Down
length(intersect(genSym.down.hk_vs_if,
                 genSym.down.ui_vs_if))
down_common <- intersect(genSym.down.hk_vs_if,
                         genSym.down.ui_vs_if)
# only HK down
length(setdiff(genSym.down.hk_vs_if,
               genSym.down.ui_vs_if))
# UI down
length(setdiff(genSym.down.ui_vs_if,
               genSym.down.hk_vs_if))

sets_down <- c(
  rep('HK_vs_Live',length(genSym.down.hk_vs_if)),
  rep('Uninfected_vs_Live',length(genSym.down.ui_vs_if))
)
m_ <- data.frame(
  elements = elements_down,
  sets= sets_down
)
#############
library('venneuler')
v_down <- venneuler(m_)
pdf('./figures/Final-C-VennDiagram-Down.pdf',
    width = 3,
    height = 3)
plot(v_down)
dev.off()
elements_up <- c(
  genSym.up.hk_vs_if,
  genSym.up.ui_vs_if
)
# up
length(intersect(genSym.up.hk_vs_if,
                 genSym.up.ui_vs_if))
up_common <- intersect(genSym.up.hk_vs_if,
                       genSym.up.ui_vs_if)
# HK up
length(setdiff(genSym.up.hk_vs_if,
               genSym.up.ui_vs_if))
# UI down
length(setdiff(genSym.up.ui_vs_if,
               genSym.up.hk_vs_if))

sets_up <- c(
  rep('HK_vs_Live',length(genSym.up.hk_vs_if)),
  rep('Uninfected_vs_Live',length(genSym.up.ui_vs_if))
)
#########
#########

m_up <- data.frame(
  elements = elements_up,
  sets= sets_up
)
v_up <- venneuler(m_up)
pdf('./figures/Final-C-VennDiagram-UP.pdf',
    width = 3,
    height = 3)
plot(v_up)
dev.off()

geneSym_heatmap <- c(down_common,up_common)

geneSym2geneAcc <- de_hk_vs_if_paired$geneAcc
names(geneSym2geneAcc) <- de_hk_vs_if_paired$geneSym
geneAccNgeneSym_heatmap <- geneSym2geneAcc[geneSym_heatmap]
save(geneAccNgeneSym_heatmap,
     file = './data/geneAccNgeneSym_heatmap.RData')