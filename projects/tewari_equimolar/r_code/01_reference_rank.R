library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

load("data/data_gff.rda")


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
