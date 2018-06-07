setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
library(edgeR)
source("r_code/functions.R")
theme_set(theme_bw(base_size = 9))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/bcbio/mirtop/expression_counts.tsv")

dds = DGEList(complete[, 13:ncol(complete)])
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = TRUE)

pilot = cbind(complete[, 1], counts) %>% 
    gather(sample, value, -UID) %>%
    mutate(sample = gsub("-mirbase-ready", "", sample)) %>% 
    filter(sample  %in%  meta_pilot[["fixed_name"]]) %>% 
    left_join(complete[,1:12]) %>% 
    filter(value > 0) %>% 
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>% 
    left_join(meta_pilot, by = c("sample" = "fixed_name"))


lapply(2:3, function(x){
    filter(pilot, value >= x) %>% summarize_isomir %>% 
        mutate(min_counts = x)
}) %>% bind_rows() %>% plot_summarize_isomir +
    ggsave("figures/replicates/bcbio.pdf", width = 9, height = 9)


pilot %>% expression_isomirs_by_lab_protocol_isomir %>%
    ggplot(aes(x=lab,y=counts,fill=as.factor(reps))) +
    geom_boxplot() + scale_y_log10() +
    facet_grid(lib_method_simple~isomir_type) + 
    ggsave("figures/replicates/bcbio_counts_per_isomir_type.pdf",
           width = 9, height = 9)


pilot %>% filter(!is.na(lib_method_simple), lab != "Lab1",
                 Variant == "Reference", value > 0) %>% 
    group_by(miRNA, lab, replicate, lib_method_simple) %>%
    summarise(counts = sum(value)) %>%
    group_by(miRNA, lab, lib_method_simple) %>% 
    summarise(reps = length(replicate), counts = sum(counts)) %>% 
    group_by(lab, reps, lib_method_simple) %>%
    summarise(n_isomirs = n(), counts = sum(counts)) %>%
    group_by(lab, lib_method_simple) %>%
    arrange(lib_method_simple, lab, desc(reps)) %>%
    mutate(n_mirs_cum = cumsum(n_isomirs)/sum(n_isomirs),
           counts_cum = cumsum(counts)/sum(counts)) %>% 
    ggplot(aes(color=as.factor(reps), x=n_mirs_cum, y=counts_cum,
               shape=as.factor(lab))) +
    geom_point() +
    scale_color_brewer("common:n_replicates", palette = "Set2") +
    scale_shape_discrete("laboratory") +
    facet_wrap(~lib_method_simple) + 
    xlab("% of sequences detected compared to a single replicate") +
    ylab("% of counts detected compared to a single replicated;")


pilot %>% plot_isoadd_position_by_protocol_by_lab(., "iso_add") +
    ggsave("figures/replicates/bcbio_isoadd_by_protocol_by_lab.pdf")

pilot %>% plot_isoadd_position_by_protocol_by_lab(., "iso_5p") +
    ggsave("figures/replicates/bcbio_iso5p_by_protocol_by_lab.pdf")

pilot %>% plot_isoadd_position_by_protocol_by_lab(., "iso_3p") +
    ggsave("figures/replicates/bcbio_iso3p_by_protocol_by_lab.pdf")

iso = "iso_5p"
filter(pilot, !is.na(lib_method_simple), lab != "Lab1") %>% 
    filter(!!sym(iso) != 0) %>% 
    select(!!sym(iso),
           miRNA, replicate, lab, value, replicate, lib_method_simple) %>% 
    gather(isomir_type, size,
           -value, -miRNA, -lab, -replicate, -lib_method_simple ) %>% 
    filter(size != 0) %>% 
    group_by(miRNA, lab, replicate, lib_method_simple, size) %>%
    summarise(counts = sum(value)) %>%
    group_by(miRNA, lab, lib_method_simple, size) %>% 
    summarise(reps = length(replicate), counts = sum(counts)) %>% 
    group_by(lab, reps, lib_method_simple, size) %>%
    summarise(n_isomirs = n(), counts = sum(counts)) %>%
    group_by(lab, lib_method_simple, size) %>%
    arrange(size, lib_method_simple, lab, desc(reps)) %>%
    mutate(n_isomirs_cum = cumsum(n_isomirs),
           counts_cum = cumsum(counts)) %>% 
    ggplot(aes(x = reps, y = n_isomirs_cum, group = lab, color = lab)) +
    geom_point() +
    geom_line() +
    facet_grid(size~lib_method_simple)
