library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

prepare =  . %>% filter(ref_is_1 == 1) %>%
    dplyr::count(protocol, sample, pct_cat, iso) %>% 
    group_by(protocol) %>% 
    mutate(sample2 = paste0(protocol, "_", 1:length(unique(sample)))) %>% 
    group_by(sample2, protocol) %>% 
    mutate(pct_total = n/sum(n)*100,
           iso = ifelse(grepl("snp ", iso), "snp + other", iso),
           iso = ifelse(grepl("add3p ", iso), "add3p + other", iso))

bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "tewari synthetic"),
    
    vandijk  %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "vandijk"), 
    
    custom  %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "tewari_custom") 
) %>% 
    ggplot(aes(x = protocol, y = pct_total, color = pct_cat)) +
    geom_boxplot(outlier.color = NA) +
    facet_grid(iso~tool, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("PCT") +
    scale_color_manual("IMPORTANCE",
                       values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggsave("results/04_other_data/04_other_data_pct_g1.pdf", height = 7)


bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "tewari"),
    
    vandijk  %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "vandijk"), 
    
    custom  %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "tewari_custom"),
    
    dsrg  %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        ungroup() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "dsrg",
               protocol = ifelse(protocol == "ilmn", "tru", protocol)),
    plasma  %>% 
        filter(pct > 1) %>% 
        prepare() %>%
        ungroup() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "human plasma",
               protocol = ifelse(protocol == "TrueSeq", "tru", "neb"))
) %>% 
    ggplot(aes(x = protocol, y = pct_total, color = pct_cat)) +
    geom_boxplot(outlier.color = NA) +
    facet_grid(iso~tool, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("PCT") +
    scale_color_manual("IMPORTANCE",
                       values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggsave("results/04_other_data/04_other_alldata_pct_g1.pdf", height = 9)

bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "tewari"),
    plasma  %>% 
        filter(pct > 1) %>% 
        prepare() %>%
        ungroup() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "human plasma",
               protocol = ifelse(protocol == "TrueSeq", "tru", "neb"))
) %>% 
    ggplot(aes(x = protocol, y = pct_total, color = pct_cat)) +
    geom_boxplot(outlier.color = NA) +
    facet_grid(iso~tool, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("PCT") +
    scale_color_manual("IMPORTANCE",
                       values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggsave("results/04_other_data/04_other_alltewari_pct_g1.pdf", height = 9)


bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "tewari"),
    plasma  %>% 
        filter(pct > 1) %>% 
        prepare() %>%
        ungroup() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "human plasma",
               protocol = ifelse(protocol == "TrueSeq", "tru", "neb"))
) %>% filter(protocol=="tru") %>% 
    ggplot(aes(pct_total, color = tool)) +
    geom_density(stat = "ecdf") +
    facet_grid(pct_cat~iso, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("PCT") # +
    #scale_color_manual("IMPORTANCE",
    #                   values = RColorBrewer::brewer.pal(7, "Dark2")[3:7])