plot_stats = function(df){
    ggplot(df, aes(x = lab, y = counts, fill = as.factor(replicate))) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_wrap(~category, scales = "free_y", nrow=4) +
        scale_fill_brewer("replicates", palette = "Set2")
}

save_stats = function(stats, tool){

    stats %>% filter(grepl("_count", category),
                     lib_method_simple == "TruSeq") %>% 
        plot_stats() +
        ggtitle(paste0(tool, " - TrueSeq")) +
        ggsave(paste0("figures/stats/", tool, "_truseq_count.png"))
    stats %>% filter(grepl("_count", category),
                     lib_method_simple == "NEBNext") %>% 
        plot_stats() +
        ggtitle(paste0(tool, " - NEBNext")) +
        ggsave(paste0("figures/stats/", tool, "_nebnext_count.png"))
    
    
    stats %>% filter(grepl("_sum", category),
                     lib_method_simple == "TruSeq") %>% 
        plot_stats() +
        ggtitle(paste0(tool, " - TrueSeq")) +
        ggsave(paste0("figures/stats/", tool, "_trueseq_sum.png"))
    stats %>% filter(grepl("_sum", category),
                     lib_method_simple == "NEBNext") %>% 
        plot_stats() +
        ggtitle(paste0(tool, " - NEBNext")) +
        ggsave(paste0("figures/stats/", tool, "_nebnext_sum.png"))
    
}

summarize_isomir = . %>% filter(!is.na(lib_method_simple), lab != "Lab1") %>% 
    filter(iso_5p != 0 | iso_3p != 0 | iso_add_nt != 0 | iso_snp_nt != 0) %>% 
    select(iso_5p, iso_3p, iso_add_nt, iso_snp_nt,
           miRNA, replicate, lab, value, replicate, lib_method_simple) %>% 
    gather(isomir_type, size,
           -value, -miRNA, -lab, -replicate, -lib_method_simple ) %>%
    filter(size != 0) %>% 
    group_by(miRNA, isomir_type, lab, replicate, lib_method_simple, size) %>%
    summarise(counts = sum(value)) %>%
    group_by(isomir_type, miRNA, lab, lib_method_simple, size) %>% 
    summarise(reps = length(replicate), counts = sum(counts)) %>% 
    group_by(lab, reps, lib_method_simple, isomir_type) %>%
    summarise(n_isomirs = n(), counts = sum(counts)) %>%
    group_by(lab, lib_method_simple, isomir_type) %>%
    arrange(isomir_type, lib_method_simple, lab, desc(reps)) %>%
    mutate(n_isomirs_cum = cumsum(n_isomirs)/sum(n_isomirs),
           counts_cum = cumsum(counts)/sum(counts))

plot_summarize_isomir = function(df) {
    ggplot(df, aes(color=as.factor(reps), x=n_isomirs_cum, y=counts_cum,
               shape=as.factor(lab),
               size=as.factor(min_counts))) +
        geom_point() +
        scale_color_brewer("common:n_replicates", palette = "Set2") +
        scale_size_discrete("filter:min_counts", range = c(1, 2.5)) +
        scale_shape_discrete("laboratory") +
        facet_grid(lib_method_simple~isomir_type) + 
        xlab("% of sequences detected compared to single sample") +
        ylab("% of counts detected compared to single sample")
}

expression_isomirs_by_lab_protocol_isomir = . %>% 
    filter(!is.na(lib_method_simple), lab != "Lab1") %>% 
    filter(iso_5p != 0 | iso_3p != 0 | iso_add_nt != 0 | iso_snp_nt != 0) %>% 
    select(iso_5p, iso_3p, iso_add_nt, iso_snp_nt,
           miRNA, replicate, lab, value, replicate, lib_method_simple) %>% 
    gather(isomir_type, size,
           -value, -miRNA, -lab, -replicate, -lib_method_simple ) %>% 
    filter(size != 0) %>% 
    group_by(miRNA, isomir_type, lab, replicate, lib_method_simple, size) %>%
    summarise(counts = sum(value)) %>%
    group_by(isomir_type, miRNA, lab, lib_method_simple, size) %>% 
    summarise(reps = length(replicate), counts = sum(counts))


plot_isoadd_position_by_protocol_by_lab = function(df, iso = NULL){
    filter(df, !is.na(lib_method_simple), lab != "Lab1") %>% 
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
        mutate(n_isomirs_cum = cumsum(n_isomirs)/sum(n_isomirs),
               counts_cum = cumsum(counts)/sum(counts)) %>% 
        ggplot(aes(color=as.factor(reps), x=n_isomirs_cum, y=counts_cum,
                   shape=as.factor(lab))) +
        geom_point() +
        scale_color_brewer("common:n_replicates", palette = "Set2") +
        scale_shape_discrete("laboratory") +
        facet_grid(lib_method_simple~size) + 
        xlab("% of sequences detected compared to single sample") +
        ylab("% of counts detected compared to single sample")
}