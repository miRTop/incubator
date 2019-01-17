setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
source("r_code/functions.R")

theme_set(theme_bw(base_size = 11))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/isomirsea/mirtop/expression_counts.tsv.gz")

dds = DGEList(complete[, 13:ncol(complete)])
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = TRUE)

pilot = cbind(complete[, 1], counts) %>%
    gather(sample, value, -UID) %>%
    mutate(sample = gsub("_tagMir-all", "", sample)) %>%
    filter(sample  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(complete[,1:12]) %>%
    filter(value >= 1) %>%
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>%
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>%
    filter(abs(iso_5p)<4, abs(iso_3p)<4, abs(iso_add)<4 )


dds$samples %>% rownames_to_column("sample") %>%
    mutate(sample = gsub("_tagMir-all", "", sample)) %>%
    filter(sample  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>%
    ggplot(aes(x = replicate, y = lib.size)) +
    geom_bar(stat = "identity") +
    facet_grid(lab~lib_method_simple) +
    ggsave("figures/replicates/isomirsea_libsize.png", width = 7, height = 9)

lapply(2:3, function(x){
    filter(pilot, value >= x) %>%
        summarize_isomir %>%
        mutate(min_counts = x)
}) %>% bind_rows() %>%
    plot_summarize_isomir +
    ggsave("figures/replicates/isomirsea.png", width = 9, height = 9)


pilot %>% expression_isomirs_by_lab_protocol_isomir %>%
    ggplot(aes(x=lab,y=counts,fill=as.factor(reps))) +
    geom_boxplot() + scale_y_log10() +
    facet_grid(lib_method_simple~isomir_type) +
    ggsave("figures/replicates/isomirsea_counts_per_isomir_type.png",
           width = 9, height = 9)
