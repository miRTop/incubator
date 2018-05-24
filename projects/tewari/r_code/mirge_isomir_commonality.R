setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_bw(base_size = 9))

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/mirge/mirtop/expression_counts.tsv")
keys = read_csv("tools/mirge/sample_fn_key.txt", col_names = c("fn", "inside")) %>% 
    mutate(fn = gsub("_isomiRs.gff", "", fn))

pilot = complete[, c(1, 13:ncol(complete))] %>%
    gather(sample, value, -UID) %>%
    inner_join(keys, by = c("sample" = "inside")) %>% 
    filter(fn  %in%  meta_pilot[["fixed_name"]]) %>% 
    left_join(complete[,1:12]) %>% 
    filter(value > 0) %>% 
    mutate(Variant = ifelse(is.na(Variant), "Reference", Variant)) %>% 
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>% 
    filter(abs(iso_5p)<4, abs(iso_3p)<4, abs(iso_add)<4 )


count_iso = . %>% filter(!is.na(lib_method_simple), lab != "Lab5") %>% 
    filter(iso_5p != 0 | iso_3p != 0 | iso_add_nt != 0 | iso_snp_nt != 0) %>% 
    select(iso_5p, iso_3p, iso_add_nt, iso_snp_nt,
           miRNA, replicate, lab, value, replicate, lib_method_simple) %>% 
    gather(isomir_type, size,
           -value, -miRNA, -lab, -replicate, -lib_method_simple ) %>% 
    group_by(miRNA, isomir_type, lab, replicate, lib_method_simple, size) %>%
    summarise(counts = sum(value)) %>%
    group_by(isomir_type, miRNA, lab, lib_method_simple, size) %>% 
    summarise(reps = length(replicate), counts = sum(counts)) %>% 
    group_by(lab, reps, lib_method_simple, isomir_type) %>%
    summarise(n_isomirs = n(), counts = sum(counts)) %>%
    group_by(lab, lib_method_simple, isomir_type) %>%
    arrange(isomir_type, lib_method_simple, lab, desc(reps)) %>%
    mutate(n_isomirs_cum = cumsum(n_isomirs), counts_cum = cumsum(counts))

lapply(1:3, function(x){
    filter(pilot, value >= x) %>% count_iso %>% 
        mutate(min_counts = x)
}) %>% bind_rows() %>% 
    ggplot(aes(color=as.factor(reps), x=n_isomirs_cum, y=counts_cum,
                   shape=as.factor(lab),
                   size=as.factor(min_counts))) +
    geom_point() +
    scale_color_brewer("common:n_replicates", palette = "Set2") +
    scale_size_discrete("filter:min_counts", range = c(1, 2.5)) +
    scale_shape_discrete("laboratory") +
    scale_y_log10() + 
    scale_x_log10() +
    facet_grid(lib_method_simple~isomir_type) + 
    ggsave("figures/replicates/mirge.pdf", width = 9, height = 9)


pilot %>% filter(!is.na(lib_method_simple), lab != "Lab5") %>% 
    filter(iso_5p != 0 | iso_3p != 0 | iso_add_nt != 0 | iso_snp_nt != 0) %>% 
    select(iso_5p, iso_3p, iso_add_nt, iso_snp_nt,
           miRNA, replicate, lab, value, replicate, lib_method_simple) %>% 
    gather(isomir_type, size,
           -value, -miRNA, -lab, -replicate, -lib_method_simple ) %>% 
    group_by(miRNA, isomir_type, lab, replicate, lib_method_simple, size) %>%
    summarise(counts = sum(value)) %>%
    group_by(isomir_type, miRNA, lab, lib_method_simple, size) %>% 
    summarise(reps = length(replicate), counts = sum(counts)) %>%
    ggplot(aes(x=lab,y=counts,fill=as.factor(reps))) +
    geom_boxplot() + scale_y_log10() +
    facet_grid(lib_method_simple~isomir_type) + 
    ggsave("figures/replicates/mirge_counts_per_isomir_type.pdf",
           width = 9, height = 9)

