setwd(here::here())
library(tidyverse)
theme_set(theme_bw())

meta_pilot = read_csv("meta_pilot.csv")

bcbio_stats = read_csv("tools/bcbio/stats/mirtop_stats.txt", skip = 1) %>% 
    mutate(sample=gsub("-mirbase-ready", "", sample)) %>% 
    inner_join(meta_pilot,
               by = c("sample" = "fixed_name"))

library(ggplot2) 
ggplot(bcbio_stats %>% filter(grepl("_count", category),
                              lib_method_simple == "TruSeq"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("bcbio - TrueSeq") +
    ggsave("figures/stats/bcbio_truseq_count.png")
ggplot(bcbio_stats %>% filter(grepl("_count", category),
                              lib_method_simple == "NEBNext"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("bcbio - NEBNext") +
    ggsave("figures/stats/bcbio_nebnext_count.png")


ggplot(bcbio_stats %>% filter(grepl("_sum", category),
                              lib_method_simple == "TruSeq"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("bcbio - TrueSeq") +
    ggsave("figures/stats/bcbio_trueseq_sum.png")
ggplot(bcbio_stats %>% filter(grepl("_sum", category),
                              lib_method_simple == "NEBNext"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("bcbio - NEBNext") +
    ggsave("figures/stats/bcbio_nebnext_sum.png")
