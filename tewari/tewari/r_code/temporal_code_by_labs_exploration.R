
add_truseq = pilot %>% filter_labs %>% 
    filter(!!sym("iso_add_nt") != 0, lib_method_simple=="TrueSeq") %>%
    select(UID, Variant,
           miRNA, sample, replicate, lab, value, replicate, lib_method_simple) %>%
    group_by(UID, sample) %>% 
    distinct(.keep_all = TRUE) %>%
    group_by(UID, lab) %>%
    mutate(n_replicate=length(unique(replicate))) %>%
    group_by(UID, n_replicate) %>% 
    mutate(n_labs = length(unique(lab))) %>% 
    arrange(desc(n_labs))  %>%
    group_by(sample) %>% 
    arrange(sample, desc(value)) %>% 
    mutate(pct=value/sum(value),
           rank = 1:n())
    

add_truseq  %>% group_by(n_labs, n_replicate) %>% 
    summarise(counts=sum(value),
              sequences=n()) %>% 
    mutate(total = sum(sequences),
           pct = sequences/total)

# boxplot to show expression of different overlap
add_truseq %>%
    ggplot(aes(x=as.factor(n_labs), y=log2(value), color=as.factor(n_replicate))) +
    geom_boxplot() +
    ylim(0,20)

add_truseq %>% 
    mutate(t3=ifelse(grepl("iso_3p", Variant),"3p","0"),
           t5=ifelse(grepl("iso_5p", Variant),"5p","0"),
           snp=ifelse(grepl("snp", Variant),"snp","0")) %>%
    unite("simple_variant",t3, t5, snp, sep="_") %>% 
    ggplot(aes(x=simple_variant)) +
    geom_bar() +
    facet_grid(n_labs~n_replicate) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# scatter plot to show expression and pct of importance of different overlap
add_truseq %>%
    ggplot(aes(x=pct, y=log2(value))) +
    geom_point() +
    facet_grid(n_replicate~n_labs)  

# How many we detect from the most accurate group only using counts
add_truseq  %>% filter(value > 2^5) %>%
    group_by(n_labs, n_replicate) %>% 
    summarise(counts=sum(value),
              sequences=n()) %>% 
    ungroup() %>% 
    mutate(total = sum(sequences),
           pct = sequences/total)

# If we use top 5000 sequences, would that recover that group?
add_truseq %>%  filter(rank<100) %>% 
    ungroup %>% 
     dplyr::count(n_labs,n_replicate) %>% 
    mutate(total = sum(n),
           pct = n/total)
add_truseq %>%  filter(rank<500) %>% 
    ungroup %>% 
    dplyr::count(n_labs,n_replicate) %>% 
    mutate(total = sum(n),
           pct = n/total)
add_truseq %>%  filter(rank<1000) %>% 
    ungroup %>% 
    dplyr::count(n_labs,n_replicate) %>% 
    mutate(total = sum(n),
           pct = n/total)


# which sequences are those one
add_truseq %>% filter(n_labs==3, n_replicate==4) %>% 
    ungroup() %>% 
    select(UID, sample, value) %>%
    distinct() %>% 
    spread(sample, value) %>% head %>% as.data.frame

