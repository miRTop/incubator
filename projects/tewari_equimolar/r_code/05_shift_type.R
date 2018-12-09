library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

prepare =  . %>% filter(ref_is_1 == 1) %>%
    dplyr::count(protocol, sample, pct_cat, iso_shift_nt) %>% 
    group_by(protocol) %>% 
    mutate(sample2 = paste0(protocol, "_", 1:length(unique(sample)))) %>% 
    group_by(sample2, protocol) %>% 
    mutate(pct_shift = n/sum(n)*100)


equimolar_razer3 %>% 
    filter(ref_is_1 == 1, iso_shift_nt!="0", (iso=="shift5p" | iso=="shift3p")) %>%
    prepare() %>% 
    ggplot(aes(x = protocol, y = pct_shift, color = pct_cat)) +
    geom_boxplot(fill = NA) +
    facet_wrap(~iso_shift_nt, nrow = 4) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggtitle("shift type") +
    ylab("% of isomiRs with shifting events") +
    scale_color_manual(values = RColorBrewer::brewer.pal(7, "Dark2")) +
    ggsave("results/05_shift_type/05_shift_type.pdf", height = 7)


equimolar_razer3 %>% 
    filter(ref_is_1 == 1) %>% 
    group_by(sample) %>% 
    mutate(n_mi_rna = length(unique(mi_rna))) %>% 
    filter(iso_shift_nt!="0", (iso=="shift5p" | iso=="shift3p")) %>%
    group_by(sample, n_mi_rna, protocol) %>% 
    summarize(n_shift = length(unique(mi_rna))) %>% 
    distinct(sample, protocol, n_mi_rna, n_shift) %>% 
    mutate(pct = n_shift / n_mi_rna * 100) %>% 
    ggplot(aes(sample, pct, fill = protocol)) +
    geom_bar(stat="identity") +
    ylab("% of miRNAs with shift events") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    ggsave("results/05_shift_type/05_shift_mirna_pct.pdf", width = 7)


equimolar_razer3 %>% 
    filter(ref_is_1 == 1) %>% 
    group_by(sample) %>% 
    mutate(n_mi_rna = length(unique(mi_rna))) %>% 
    filter(iso_shift_nt!="0", (iso=="shift5p" | iso=="shift3p"),
           iso_shift_nt  %in% c("-1 0", "0 -1")) %>%
    group_by(sample, n_mi_rna, protocol) %>% 
    summarize(n_shift = length(unique(mi_rna))) %>% 
    distinct(sample, protocol, n_mi_rna, n_shift) %>% 
    mutate(pct = n_shift / n_mi_rna * 100) %>% 
    ggplot(aes(sample, pct, fill = protocol)) +
    geom_bar(stat="identity") +
    ylab("% of miRNAs with shift events") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    ggsave("results/05_shift_type/05_shift1_mirna_pct.pdf", width = 7)
