library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

prepare =  . %>% filter(ref_is_1 == 1) %>%
    dplyr::count(protocol, sample, short, pct_cat, iso) %>% 
    group_by(protocol) %>% 
    mutate(sample2 = paste0(protocol, "_", 1:length(unique(sample)))) %>% 
    group_by(sample2, protocol) %>% 
    mutate(pct_total = n/sum(n)*100,
           iso = gsub("snp", "snv", iso),
           iso = ifelse(grepl("snv ", iso), "snv + other", iso),
           iso = ifelse(grepl("add3p ", iso), "add3p + other", iso)) %>% 
    filter(iso  %in% c("shift5p", "shift3p", "snv", "reference"))

bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "tewari synthetic"),
    
    vandijk  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "vandijk"), 
    
    custom  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "tewari_custom"),
    
    narrykim  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "narrykim_synthetic"),
) %>% 
    ggplot(aes(x = protocol, y = pct_total, color = pct_cat)) +
    geom_boxplot(outlier.color = NA) +
    facet_grid(iso~tool, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("% of sequences") +
    scale_color_manual("IMPORTANCE",
                       values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggsave("results/04_other_data/04_other_data_pct_g5.pdf", height = 7)


bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "tewari"),
    
    vandijk  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "vandijk"), 
    
    custom  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "tewari_custom"),
    
    dsrg  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        ungroup() %>% 
        mutate(tool = "dsrg",
               protocol = ifelse(protocol == "ilmn", "tru", protocol)),
    
    narrykim  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "narrykim_synthetic"),

    plasma  %>% 
        filter(pct > 5) %>% 
        prepare() %>%
        ungroup() %>% 
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
    ggsave("results/04_other_data/04_other_alldata_pct_g5.pdf", height = 9)

# barplots
bind_rows(
    equimolar_razer3 %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "tewari"),
    
    vandijk  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "vandijk"), 
    
    custom  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "tewari_custom"),
    
    dsrg  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        ungroup() %>% 
        mutate(tool = "dsrg",
               protocol = ifelse(protocol == "ilmn", "tru", protocol)),
    
    narrykim  %>% 
        filter(pct > 5) %>% 
        prepare() %>% 
        mutate(tool = "narrykim_synthetic"),
    
    plasma  %>% 
        filter(pct > 5) %>% 
        prepare() %>%
        ungroup() %>% 
        mutate(tool = "human plasma",
               protocol = ifelse(protocol == "TrueSeq", "tru", "neb"))
) %>% 
    ungroup() %>% 
    mutate(index = as.numeric(as.factor(sample)),
           clean_name = paste0(tool, "_", protocol, "_", index)) %>% 
    ggplot(aes(x = clean_name, y = pct_total, fill = pct_cat)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(size=7, angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("% of sequences") +
    ylim(0,100) +
    scale_fill_manual("IMPORTANCE",
                       values = RColorBrewer::brewer.pal(7, "Dark2")[4:7]) +
    facet_wrap(~iso, nrow = 4) +
    ggsave("results/04_other_data/04_other_alldata_pct_g5_bar.pdf", height = 9, width = 12)

bind_rows(narrykim %>% mutate(short=paste("nk_", short)),
          narrykim_human %>% mutate(short=paste("nkh_", short)))%>%
    filter(iso  %in% c("shift5p", "shift3p", "snp")) %>%
    ggplot(aes(y=short,x=pct)) +
    ggridges::geom_density_ridges(color=NA) +
    geom_vline(xintercept = 10) +
    facet_wrap(~iso, nrow = 1,  scales = "free_x") +
    ggsave("results/04_other_data/04_other_nk_pct_g10_density.pdf", height = 9)

bind_rows(plasma %>% mutate(short=paste("twh_", short)),
          equimolar_razer3 %>% mutate(short=paste("tw_", short)))%>%
    filter(iso  %in% c("shift5p", "shift3p", "snp")) %>%
    ggplot(aes(y=short,x=pct)) +
    ggridges::geom_density_ridges(color=NA) +
    geom_vline(xintercept = 10) +
    facet_wrap(~iso, nrow = 1,  scales = "free_x") +
    ggsave("results/04_other_data/04_other_tw_pct_g10_density.pdf", height = 9)



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