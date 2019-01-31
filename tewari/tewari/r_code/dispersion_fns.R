setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
source("r_code/functions.R")
theme_set(theme_bw(base_size = 14))


get_dispersion <- function(c){
    c = c[rowSums(c)>0,]
    
    d = DGEList(c) %>%
        calcNormFactors(.) %>%
        estimateGLMCommonDisp(.) %>% 
        estimateGLMTagwiseDisp(.) 
    
    data.frame(feature = rownames(c),
               avg = d[["AveLogCPM"]],
               dispersion = d[["tagwise.dispersion"]],
               stringsAsFactors = F)
}

analysis = . %>% 
    group_by(miRNA, group) %>% 
    mutate(total = sum(avg), pct = avg / total *100) %>% 
    mutate(pct_cat = cut(pct,
                         breaks = c(-1, 0.1, 1, 5, 10, 20, 50, 101),
                         labels = c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")))

annotation = . %>% 
    mutate(fix_add = iso_add ==  -1 * iso_3p & iso_add != 0,
           iso_add = ifelse(fix_add, 0, iso_add),
           iso_3p = ifelse(fix_add, 0, iso_3p),
           iso_snp = ifelse(fix_add, 1, iso_snp)) %>% 
    mutate(is_snp = ifelse(iso_snp != 0, "snp", ""),
           is_add = ifelse(iso_add != 0, "add3p", ""),
           is_5p = ifelse(iso_5p != 0, "shift5p", ""),
           is_3p = ifelse(iso_3p != 0, "shift3p", ""),
    ) %>% 
    mutate(is_loss_5p = ifelse(iso_5p < 0, "loss_5p", ""),
           is_loss_5p = ifelse(iso_5p > 0, "gain_5p", is_loss_5p),
           is_loss_3p = ifelse(iso_3p < 0, "loss_3p", ""),
           is_loss_3p = ifelse(iso_3p > 0, "gain_3p", is_loss_3p)) %>%
    unite("iso_shift_nt", iso_5p, iso_3p, sep = " ", remove = FALSE) %>%
    unite("iso_nt", iso_snp, iso_add, iso_5p, iso_3p, sep = "_") %>%
    #unite("iso_detail", iso_add_nt, iso_5p_nt, iso_3p_nt, sep = "_") %>%
    unite("iso", is_snp, is_add, is_5p, is_3p, sep = " ") %>%
    mutate(iso = trimws(iso),
           iso = gsub("[ ]+", " ", iso),
           iso = ifelse(iso=="", "reference",iso)) %>% 
    unite("iso_shift", is_loss_5p, is_loss_3p, sep =  " ") %>%
    mutate(iso_shift = trimws(iso_shift),
           iso_shift = ifelse(iso_shift=="", "no-loss",iso_shift)) %>% 
    select(UID, miRNA, iso, iso_shift, iso_shift_nt) 
