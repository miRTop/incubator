setwd(here::here())

meta_pilot = read_csv("meta_pilot.csv")

## needs to get name inside the file and name of the file to map to original
keys = read_csv("tools/mirge/sample_fn_key.txt", col_names = c("fn", "inside")) %>% 
    mutate(fn = gsub("_isomiRs.gff", "", fn))
mirge_stats = read_csv("tools/mirge/stats/mirtop_stats.txt", skip = 1) %>% 
    inner_join(keys, by = c("sample" = "inside")) %>% 
    inner_join(meta_pilot,
               by = c("fn" = "fixed_name")) %>% 
    mutate(sample = fn)

library(ggplot2) 
ggplot(mirge_stats %>% filter(grepl("_count", category),
                              lib_method_simple == "TruSeq"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("miRge - TrueSeq") +
    ggsave("figures/stats/mirge_trueseq_count.png")
ggplot(mirge_stats %>% filter(grepl("_count", category),
                              lib_method_simple == "NEBNext"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("miRge - NEBNext") +
    ggsave("figures/stats/mirge_nebnext_count.png")

ggplot(mirge_stats %>% filter(grepl("_sum", category),
                              lib_method_simple == "TruSeq"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("miRge - TrueSeq") +
    ggsave("figures/stats/mirge_trueseq_sum.png")
ggplot(mirge_stats %>% filter(grepl("_sum", category),
                              lib_method_simple == "NEBNext"),
       aes(x = lab, y = counts, fill = as.factor(replicate))) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~category, scales = "free_y", nrow=4) +
    ggtitle("miRge - NEBNext") +
    ggsave("figures/stats/mirge_nebnext_sum.png")