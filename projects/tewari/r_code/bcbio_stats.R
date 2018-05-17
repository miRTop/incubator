setwd(here::here())

meta_pilot = read_csv("meta_pilot.csv")

bcbio_stats = read_csv("tewari/bcbio/stats/mirtop_stats.txt", skip = 1) %>% 
    mutate(sample=gsub("-mirbase-ready", "", sample)) %>% 
    inner_join(meta_pilot,
               by = c("sample" = "fixed_name"))

library(ggplot2) 
ggplot(bcbio_stats %>% filter(grepl("_count", category),
                              lib_method_simple == "TruSeq"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggsave("figures/stats/bcbio_count.pdf")

ggplot(bcbio_stats %>% filter(grepl("_sum", category),
                              lib_method_simple == "TruSeq"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggsave("figures/stats/bcbio_sum.pdf")
