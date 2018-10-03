library(isomiRs)
meta_pilot = read_csv("../tewari/meta_pilot.csv") %>% 
    as.data.frame()
rownames(meta_pilot) = paste0(meta_pilot$fixed_name, "-mirbase-ready.counts")
fns = file.path("../tewari/tools/bcbio", meta_pilot$fixed_name, rownames(meta_pilot))
iso_plasma = IsomirDataSeqFromFiles(fns, meta_pilot)
save(iso_plasma, file = "data/iso_plasma.rda")

library(tidyverse)
load("results/ranked.rda")

bcbio = metadata(iso_plasma)[["rawData"]]
meta_pilot = read_csv("../tewari/meta_pilot.csv")

library(edgeR)
dds = DGEList(bcbio[, 7:ncol(bcbio)])
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = FALSE)


plasma = bcbio %>%
    mutate(is_snp = ifelse(mism != 0, "snp", ""),
           is_add = ifelse(add != 0, "add", ""),
           is_t5 = ifelse(t5 != 0, "t5", ""),
           is_t3 = ifelse(t3 != 0, "t3", "")) %>% 
    mutate(is_loss_t5 = ifelse(t5 != 0 & grepl("[a-z]+", t5), "loss_t5", "NaN"),
           is_loss_t5 = ifelse(t5 != 0 & grepl("[A-Z]+", t5), "gain_t5", is_loss_t5),
           is_loss_t3 = ifelse(t3 != 0 & grepl("[a-z]+", t3), "loss_t3", "NaN"),
           is_loss_t3 = ifelse(t3 != 0 & grepl("[A-Z]+", t3), "gain_t3", is_loss_t3)) %>%
    unite("iso_nt", mism, add, t5, t3, sep = "_") %>%
    unite("iso", is_snp, is_add, is_t5, is_t3, sep = ".") %>%
    unite("iso_loss", is_loss_t5, is_loss_t3, sep = ".") %>%
    gather(sample, value, -mir, -iso, -iso_nt, -iso_loss, -seq) %>% 
    mutate(sample = gsub("-mirbase-ready.counts", "", sample)) %>%
    left_join(meta_pilot, by = c("sample" = "fixed_name")) %>%
    group_by(sample, mir) %>% 
    filter(value > 0) %>% 
    arrange(sample, mir, desc(value)) %>%
    mutate(rank = 1:n(),
           total = sum(value),
           pct = value / total * 100) %>% 
    ungroup() %>% 
    group_by(mir, iso, lab) %>% 
    mutate(reproducible_lab = n()) %>%
    group_by(mir, iso) %>% 
    mutate(reproducible_protocol = length(unique(lib_method_simple))) %>%
    group_by(mir, sample) %>%
    mutate(ref_is_1 = length(mir[rank == 1 & iso == "..."])) %>% 
    ungroup %>% 
    group_by(sample) %>%
    mutate(pct_cat = cut(pct,
                         breaks = c(-1, 0.1, 1, 5, 10, 20, 50, 101),
                         labels = c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")),
           mir_coverage = cut(total,
                              breaks = c(-1, 1, 10, 100, 1000, 1e30),
                              labels = c("<1", "1-10", "10-100", "100-1000", "<1000"))) %>% 
    ungroup() %>% 
    filter(mir  %in% ranked$mir) %>% 
    unite("short", c("lib_method_simple", "lab", "replicate"), remove = F)



full_join(
    dplyr::count(distinct(plasma, sample, mir), sample),
    dplyr::count(filter(plasma, iso == "...", rank == 1), sample),
    by = "sample", suffix = c("_total", "_ref_is_1")
) %>% mutate(pct = n_ref_is_1/n_total*100)


plasma %>% 
    filter(ref_is_1 == 1) %>%
    ggplot(aes(x = short, fill = pct_cat)) +
    geom_bar() +
    facet_wrap(~iso, nrow = 4, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plasma %>% 
    filter(ref_is_1 == 1) %>%
    ggplot(aes(x = short, fill = pct_cat)) +
    geom_bar() +
    facet_wrap(~iso, nrow = 4) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


plasma %>% 
    filter(ref_is_1 == 1, !grepl("snp", iso)) %>%
    ggplot(aes(x = short, fill = pct_cat)) +
    geom_bar() +
    facet_grid(mir_coverage~iso) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plasma %>% 
    filter(ref_is_1 == 1, !grepl("snp", iso)) %>%
    group_by(sample) %>% 
    mutate(total_isomirs = n()) %>% 
    group_by(lib_method_simple, sample, iso, pct_cat) %>% 
    summarize(pct_sequences = n()/total_isomirs[1]) %>% 
    ggplot(aes(x = iso, color = pct_cat, y=pct_sequences)) +
    geom_boxplot() + 
    facet_wrap(~lib_method_simple) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


plasma %>%
    filter(ref_is_1 != 1, rank == 1) %>% 
    ggplot(aes(short)) +
    geom_bar() +
    facet_wrap(~iso) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
