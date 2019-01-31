library(DESeq2)
library(edgeR)

ddslow = makeExampleDESeqDataSet(n = 20000, m = 4, dispMeanRel = function(x) 4/x + 0.1) %>% counts(.)
pheatmap::pheatmap(cor(ddslow), breaks = seq(0,1.1,0.1), color = RColorBrewer::brewer.pal(10, "Spectral"))
ddslow %>% DGEList(.) %>%  calcNormFactors(.) %>% estimateGLMTrendedDisp(.)  %>% .[["trended.dispersion"]] %>% median

ddshigh = makeExampleDESeqDataSet(n = 20000, m = 4, dispMeanRel = function(x) 4/x + 0.8) %>% counts(.)
pheatmap::pheatmap(cor(ddshigh), breaks = seq(0,1.1,0.1), color = RColorBrewer::brewer.pal(10, "Spectral"))
ddshigh %>% DGEList(.) %>%  calcNormFactors(.) %>% estimateGLMTrendedDisp(.)  %>% .[["trended.dispersion"]] %>% median
