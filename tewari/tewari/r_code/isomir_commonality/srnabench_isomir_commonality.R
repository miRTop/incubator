setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
source("r_code/functions.R")

theme_set(theme_bw(base_size = 11))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/srnabench/mirtop/expression_counts.tsv.gz")

# srnabench is not compatible with --add-extra option, so the samples start at column 9
dds = DGEList(complete[, 9:ncol(complete)])
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = TRUE)

pilot = cbind(complete[, 1], counts) %>%
    gather(sample, value, -UID) %>%
    filter(sample  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(complete[,1:12]) %>%
    filter(value >= 1) %>%
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>%
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>%
    filter(abs(iso_5p)<4, abs(iso_3p)<4, abs(iso_add)<4 )


dds$samples %>% rownames_to_column("sample") %>%
    filter(sample  %in%  meta_pilot[["fixed_name"]]) %>%
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>%
    ggplot(aes(x = replicate, y = lib.size)) +
    geom_bar(stat = "identity") +
    facet_grid(lab~lib_method_simple) +
    ggsave("figures/replicates/srnabench_libsize.png", width = 7, height = 9)

add_col = "iso_add"
snp_col = "iso_snp"

pilot %>%
    summarize_isomir %>%
    plot_summarize_isomir +
    ggsave("figures/replicates/srnabench.png", width = 11, height = 9)

pilot %>%
    summarize_isomirs_by_lab() %>%
    plot_summarize_isomir_by_lab() +
    ggsave("figures/labs/srnabench.png", width = 11, height = 9)


pilot %>% expression_isomirs_by_lab_protocol_isomir %>%
    ggplot(aes(x=lab,y=counts,fill=as.factor(reps))) +
    geom_boxplot() + scale_y_log10() +
    facet_grid(lib_method_simple~isomir_type) +
    ggsave("figures/replicates/srnabench_counts_per_isomir_type.png",
           width = 9, height = 9)
