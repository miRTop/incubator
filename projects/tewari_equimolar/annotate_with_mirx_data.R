library(tidyverse)
library(seqinr)
mirbase = read.fasta(file = "data/mature.fa.gz", as.string = TRUE, forceDNAtolower = FALSE)

mirbase_df = data.frame(id = names(mirbase),
                        sequence = gsub("U", "T", unlist(mirbase)),
                        stringsAsFactors = F)

mirx = rio::import("data/nbt.4183-S3.xlsx", skip=2) %>% 
    janitor::clean_names()

mirx_labeled = left_join(mirx, filter(mirbase_df, grepl("hsa-", id)),
          by = c("sequence")) %>% 
    mutate(id=ifelse(is.na(id), equimolar_sequence_id, id)) %>% 
    select(sequence, id)

# number of annotated sequence with mirbase in human
sum(grepl("hsa-", mirx_labeled$id))

library(bcbioSmallRna)
library(ggplot2)
library(edgeR)
load("data/bcb.rda")
ids = bcbio(bcb, "isomirs")

counts = metadata(ids)[["rawData"]][,7:27] %>% 
    as.matrix()
dds = DGEList(counts)
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = FALSE)


# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
clean_data = metadata(ids)[["rawData"]] %>% 
    left_join(mirx_labeled, by = c("seq" = "sequence")) %>% # 31493 rows
    distinct(seq, mir, .keep_all = TRUE) # 31443 rows, need to identify what sequences are lost here

metadata(ids)[["rawData"]][,1:6] %>% 
    right_join(mirx_labeled, by = c("seq" = "sequence")) %>% 
    filter(is.na(mir)) %>% .[["id"]]

normalized = cbind(metadata(ids)[["rawData"]][,1],
      counts) %>% 
    gather(sample, normalized,  -seq)
    

# already missing
setdiff(mirx_labeled$id,clean_data$mir) %>% .[grepl("hsa-",.[])]

# rank each sequence by each miRNA, and sample
parsed =  clean_data %>%
    filter(mir == id | is.na(id)) %>% # only allow families where annotation is equal to expected 31312 rows
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
    gather(sample, value, -mir, -iso, -iso_nt, -iso_loss, -id, -seq) %>%
    group_by(mir, sample) %>% # remove miRNAs that have more than one spike in in the family
    mutate(any_in_mirx = sum(grepl("hsa-", id))) %>% 
    ungroup() # 31312 seq,mir rows

# already missing
setdiff(mirx_labeled$id,parsed$mir) %>% .[grepl("hsa-",.[])]

ranked = parsed %>% filter(any_in_mirx == 1) %>% # 25665 (seq, mir rows), only allow families where the reference appears once
    left_join(normalized, by = c("sample", "seq")) %>% 
    mutate(lab=stringr::str_extract(sample,"Lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$")) %>% 
    distinct() %>% 
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
    mutate(reproducible_protocol = length(unique(protocol))) %>%
    group_by(mir, sample) %>%
    mutate(ref_is_1 = length(mir[rank == 1 & iso == "..."])) %>% 
    ungroup %>% 
    group_by(sample) %>%
    mutate(pct_cat = cut(pct,
                          breaks = c(-1, 0.1, 1, 5, 10, 20, 50, 101),
                          labels = c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")),
           mir_coverage = cut(total,
                              breaks = c(-1, 1, 10, 100, 1000, 1e30),
                              labels = c("<1", "1-10", "10-100", "100-1000", ">1000"))) %>% 
    ungroup()

save(ranked, file = "results/ranked.rda")

# We missed these families
setdiff(mirx_labeled$id,ranked$mir) %>% .[grepl("hsa-",.[])]
