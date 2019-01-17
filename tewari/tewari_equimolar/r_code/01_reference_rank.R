library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")


full_join(
    dplyr::count(distinct(equimolar_razer3, sample, mi_rna, protocol), sample, protocol),
    dplyr::count(filter(equimolar_razer3,  rank == 1, !is.na(id), iso == "reference"), sample),
    by = "sample", suffix = c("_total", "_ref_is_1")
) %>% mutate(pct = n_ref_is_1/n_total*100,
             ytext = pct + pct*0.1,
             text = paste0(round(n_total, digits = 0))) %>% 
    ggplot(aes(sample, pct, fill = protocol)) +
    geom_bar(stat = "identity") +
    geom_text(aes(sample, ytext, label = text), size = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    ggsave("results/01_reference_rank/01_reference_rank_bar.pdf")


filter(equimolar_razer3,  rank == 1, is.na(id), iso != "reference") %>% 
    mutate(iso = ifelse(grepl("snp ", iso), "snp + other", iso),
           iso = ifelse(grepl("add3p ", iso), "add3p + other", iso)) %>% 
    ggplot(aes(sample, fill = iso)) + 
    geom_bar() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Set2")) +
    ylab("# of miRNAs") +
    ggsave("results/01_reference_rank/01_ref_isnot_first_types.pdf", width = 9)

# among the non-reference being the top, isomir distribution
notop = filter(equimolar_razer3, !grepl("x4n_[abcd]", sample)) %>% 
    filter(  rank == 1, iso != "reference", is.na(id)) %>% 
    .[["mi_rna"]] %>% unique

df = filter(equimolar_razer3, !grepl("x4n_[abcd]", sample))  %>%  
    filter(rank == 1, any_in_mirx >0, mi_rna  %in% notop) %>% 
    mutate(iso = ifelse(grepl("snp ", iso), "snp + other", iso),
           iso = ifelse(grepl("add3p ", iso), "add3p + other", iso)) 
    
df_ma = select(df, mi_rna, iso, sample) %>% 
    mutate(iso = as.numeric(factor(iso))) %>% 
    spread(sample, iso, fill = 0) %>% 
    as.data.frame() %>% 
    column_to_rownames("mi_rna") %>% 
    as.matrix()
hr = hclust(dist(df_ma), method = "ward.D")
hc = hclust(dist(t(df_ma)), method = "ward.D")
sample_order = hc$labels[hc$order]
mirna_order = hr$labels[hr$order]

df %>% 
    mutate(sample = factor(sample, levels = sample_order),
              mi_rna = factor(mi_rna, levels = mirna_order),
              iso = relevel(factor(iso), "reference")) %>% 
    ggplot(aes(sample, mi_rna, fill = iso)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank()) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(8, "Set2"), na.value = "white")+
    ggsave("results/01_reference_rank/01_ref_isnot_first_heatmap.pdf", width = 9)


# counts  > 1000
notop1000 = filter(equimolar_razer3, !grepl("x4n_[abcd]", sample)) %>% 
    filter(rank == 1, iso != "reference", is.na(id), total > 1000) %>% 
    .[["mi_rna"]] %>% unique

df %>% filter(mi_rna  %in%  notop1000) %>% 
    mutate(sample = factor(sample, levels = sample_order),
           mi_rna = factor(mi_rna, levels = mirna_order),
           iso = relevel(factor(iso), "reference")) %>% 
    ggplot(aes(sample, mi_rna, fill = iso)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_blank()) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(8, "Set2"), na.value = "white") +
    ggsave("results/01_reference_rank/01_ref_isnot_first_1000_heatmap.pdf", width = 9)

    