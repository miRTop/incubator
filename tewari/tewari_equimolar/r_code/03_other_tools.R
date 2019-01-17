library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 14L))
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
        mutate(tool = "razer3"),
    
    equimolar_mirge  %>% 
        mutate(protocol = ifelse(grepl("NEBN", protocol), "neb", protocol),
               protocol = ifelse(grepl("Tru", protocol), "tru", protocol),
               protocol = ifelse(grepl("4N", protocol), "x4n", protocol),
               protocol = ifelse(grepl("Clean", protocol), "clean", protocol)) %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "mirGe"), 
    
    equimolar  %>% 
        filter(pct > 1) %>% 
        prepare() %>% 
        filter(iso  %in% c("shift5p", "shift3p", "snp")) %>% 
        mutate(tool = "bcbio") 
) %>% 
    ggplot(aes(x = protocol, y = pct_total, color = pct_cat)) +
    geom_boxplot(outlier.color = NA) +
    facet_grid(tool~iso) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("PCT") +
    scale_color_manual("IMPORTANCE",
                       values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggsave("results/03_other_tools/03_other_tools_pct_g1.pdf", height = 9)
