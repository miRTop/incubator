setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
source("r_code/functions.R")

theme_set(theme_bw(base_size = 11))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/mirge/mirtop/expression_counts.tsv.gz")
keys = read_csv("tools/mirge/sample_fn_key.txt", col_names = c("fn", "inside")) %>%
    mutate(fn = gsub("_isomiRs.gff", "", fn))

dds = DGEList(complete[, 13:ncol(complete)])
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = TRUE)

pilot = cbind(complete[, 1], counts) %>%
    gather(sample, value, -UID) %>%
    inner_join(keys, by = c("sample" = "inside")) %>%
    filter(fn  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(complete[,1:12]) %>%
    filter(value >= 1) %>%
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>%
    left_join(meta_pilot, by = c("fn" = "fixed_name")) %>%
    filter(abs(iso_5p)<4, abs(iso_3p)<4, abs(iso_add)<4 )


dds$samples %>% rownames_to_column("sample") %>%
    inner_join(keys, by = c("sample" = "inside")) %>%
    filter(fn  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(meta_pilot, by = c("fn" = "fixed_name")) %>%
    ggplot(aes(x = replicate, y = lib.size)) +
    geom_bar(stat = "identity") +
    facet_grid(lab~lib_method_simple) +
    ggsave("figures/replicates/mirge_libsize.png", width = 7, height = 9)

pilot %>%
    summarize_isomir %>%
    plot_summarize_isomir +
    ggsave("figures/replicates/mirge.png", width = 11, height = 9)

pilot %>%
    summarize_isomirs_by_lab() %>%
    plot_summarize_isomir_by_lab() +
    ggsave("figures/labs/mirge.png", width = 11, height = 9)

pilot %>% expression_isomirs_by_lab_protocol_isomir %>%
    ggplot(aes(x=lab,y=counts,fill=as.factor(reps))) +
    geom_boxplot() + scale_y_log10() +
    facet_grid(lib_method_simple~isomir_type) +
    ggsave("figures/replicates/mirge_counts_per_isomir_type.png",
           width = 9, height = 9)
