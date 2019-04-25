library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 16L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")


prepare =  . %>% filter(ref_is_1 == 1) %>%
    dplyr::count(sample, pct_cat, iso, protocol) %>% 
    group_by(sample, protocol) %>% 
    mutate(pct_total = n/sum(n)*100,
           iso = ifelse(grepl("snp ", iso), "snp + other", iso),
           iso = ifelse(grepl("add3p ", iso), "add3p + other", iso),
           iso = gsub("shift", "iso_", iso),
           iso = gsub("snp", "snv", iso),
           iso=relevel(as.factor(iso), "reference"))

equimolar_razer3 %>% 
    prepare() %>% 
    ggplot(aes(x = sample, fill = pct_cat, y = pct_total)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    ggtitle("isomiRs importance by type") +
    ylab("% of sequences") +
    ggsave("results/02_importance/02_importance.pdf", height = 9)



equimolar_razer3 %>% 
    filter(total > 1000) %>% 
    prepare() %>% 
    ggplot(aes(x = sample, fill = pct_cat, y = pct_total)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    ggtitle("isomiRs importance by type (miRNAs > 1000 counts)") +
    ylab("% of sequences") +
    ggsave("results/02_importance/02_importance_1000.pdf", height = 9)

full_join(
equimolar_razer3 %>% 
    filter(pct > 1) %>% 
    group_by(sample, protocol) %>% 
    summarise(n=n()),
equimolar_razer3 %>% 
    filter(pct > 0) %>% 
    group_by(sample) %>% 
    summarise(n=n()),
by = "sample", suffix = c("_filtered", "_total")
) %>% 
    mutate(pct = (1-n_filtered/n_total)*100,
           ytext = pct + pct*0.1,
           text = paste0(round(n_filtered, digits = 0))) %>% 
    ggplot(aes(sample, pct, fill = protocol)) +
    geom_bar(stat = "identity") +
    geom_text(aes(sample, ytext, label = text), size = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("% removed") + 
    ggsave("results/02_importance/02_pct_filtered.pdf")


equimolar_razer3 %>% 
    filter(pct > 1) %>% 
    prepare() %>% 
    ggplot(aes(x = sample, fill = pct_cat, y = pct_total)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ylab("% of sequences") +
    ggtitle("isomiRs importance by type with pct > 1") +
    ggsave("results/02_importance/02_importance_pct_g1.pdf", height = 9)

equimolar_razer3 %>% 
    filter(total > 1000, pct > 1) %>% 
    prepare() %>% 
    ggplot(aes(x = sample, fill = pct_cat, y = pct_total)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggtitle("isomiRs importance by type (miRNAs > 1000 counts)") +
    ylab("% of sequences") +
    ggsave("results/02_importance/02_importance_1000_pct_g1.pdf", height = 9)



equimolar_razer3 %>% 
    filter(pct > 1) %>% 
    prepare() %>% 
    ggplot(aes(x = protocol, y = pct_total, color = pct_cat)) +
    geom_boxplot(outlier.color = NA) +
    facet_wrap(~iso, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("PCT") +
    scale_color_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")[3:7]) +
    ggsave("results/02_importance/02_importance_pct_g1_boxplot.pdf", height = 9)


equimolar_razer3 %>% 
    filter(pct > 5) %>% 
    prepare() %>% 
    ggplot(aes(x = sample, fill = pct_cat, y = pct_total)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2")[4:7]) +
    ggtitle("isomiRs importance by type with pct > 5") +
    ggsave("results/02_importance/02_importance_pct_g5.pdf", height = 9)


equimolar_razer3 %>% 
    filter(pct > 1, ref_is_1 == 1) %>% 
    mutate(iso = ifelse(grepl("snp ", iso), "snp + other", iso),
           iso = ifelse(grepl("add3p ", iso), "add3p + other", iso)) %>% 
    ggplot(aes(x = protocol, fill = pct_cat, y = normalized)) +
    geom_violin() +
    scale_y_log10() +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(7, "Dark2")[2:7]) +
    ggtitle("isomiRs abundance by type with pct > 1") +
    ggsave("results/02_importance/02_abundance_pct_g1.pdf", height = 9)
