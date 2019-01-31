source("r_code/dispersion_fns.R")

meta_pilot = read_csv("meta_pilot.csv")

complete = read_tsv("tools/srnabench/mirtop/expression_counts.tsv.gz")
dds = DGEList(complete[, 13:ncol(complete)])
rownames(dds) = complete$UID

meta_pilot = meta_pilot %>% 
    filter(fixed_name  %in% names(complete)) %>% 
    mutate(group = paste(lab, lib_method_simple)) 
meta_pilot[["group"]] %>% unique %>% 
    lapply(., function(l){
        d = get_dispersion(dds$counts[,meta_pilot$fixed_name[meta_pilot$group==l]])
        d[["group"]] = l
        d
    }) %>% bind_rows() %>% 
    mutate(tool = "srnabench") %>% 
    left_join(complete %>% 
                  annotation(),
              by = c("feature" = "UID")) %>% 
    analysis() -> full

dds$samples %>% as.data.frame() %>% 
    rownames_to_column("sample") %>% 
    inner_join(meta_pilot, by = c("sample" = "fixed_name")) %>% 
    select(group.y, lib.size) %>% 
    write_csv("tables/dispersion/srnabench.lib.size.csv")

full %>% 
    mutate(iso = ifelse(grepl("snp ", iso), "nv + other", iso)) %>% 
    group_by(group, iso, pct_cat, tool) %>% 
    summarise(disp = median(dispersion), exp = median(avg), n = n()) %>% 
    write_csv("tables/dispersion/srnabench.csv")
