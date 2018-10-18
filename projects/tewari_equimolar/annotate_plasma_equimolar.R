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

library(ggplot2)
library(edgeR)
load("data/bcb.rda")

# this code parse the columns to get the information as needed to plot
annotate = . %>% 
    mutate(is_snp = ifelse(iso_snp != 0, "snp", ""),
           is_add = ifelse(iso_add != 0, "add", ""),
           is_5p = ifelse(iso_5p != 0, "t5", ""),
           is_3p = ifelse(iso_3p != 0, "t3", "")) %>% 
    mutate(is_loss_5p = ifelse(iso_5p != 0 & grepl("[a-z]+", iso_5p), "loss_t5", "NaN"),
           is_loss_5p = ifelse(iso_5p != 0 & grepl("[A-Z]+", iso_5p), "gain_t5", is_loss_5p),
           is_loss_3p = ifelse(iso_3p != 0 & grepl("[a-z]+", iso_3p), "loss_t3", "NaN"),
           is_loss_3p = ifelse(iso_3p != 0 & grepl("[A-Z]+", iso_3p), "gain_t3", is_loss_3p)) %>%
    unite("iso_nt", iso_snp, iso_add, iso_5p, iso_3p, sep = "_") %>%
    unite("iso", is_snp, is_add, is_5p, is_3p, sep = ".") %>%
    unite("iso_loss", is_loss_5p, is_loss_3p, sep = ".") %>%
    select(-uid, -variant, -iso_5p_nt:-iso_snp_nt ) %>% 
    gather(sample, value, -mi_rna, -iso, -iso_nt, -iso_loss, -id, -read) %>% 
    ungroup

# quantify the times sequences are detected and % of importance
analysis = . %>% 
    group_by(sample, mi_rna) %>% 
    filter(value > 0) %>% 
    arrange(sample, mi_rna, desc(value)) %>%
    mutate(rank = 1:n(),
           total = sum(value),
           pct = value / total * 100) %>% 
    ungroup() %>% 
    group_by(mi_rna, iso, lab) %>% 
    mutate(reproducible_lab = n()) %>%
    group_by(mi_rna, iso) %>% 
    mutate(reproducible_protocol = length(unique(protocol))) %>%
    group_by(mi_rna, sample) %>%
    mutate(ref_is_1 = length(mi_rna[rank == 1 & iso == "..."])) %>% 
    ungroup %>% 
    mutate(pct_cat = cut(pct,
                          breaks = c(-1, 0.1, 1, 5, 10, 20, 50, 101),
                          labels = c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")),
           mir_coverage = cut(total,
                              breaks = c(-1, 1, 10, 100, 1000, 1e30),
                              labels = c("<1", "1-10", "10-100", "100-1000", ">1000"))) %>% 
    ungroup()


# equimolar
gff = read_tsv("tools/bcbio/mirtop/expression_counts.tsv.gz") %>% 
    janitor::clean_names()

counts = gff[,13:33] %>% 
    as.matrix()
dds = DGEList(counts)
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = FALSE)
normalized = cbind(gff[,"read"],
                   counts) %>% 
    gather(sample, normalized,  -read)

# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
clean_data = gff %>% 
    left_join(mirx_labeled, by = c("read" = "sequence")) %>% # 31493 rows
    distinct(read, mi_rna, .keep_all = TRUE) # 31443 rows, need to identify what sequences are lost here

gff %>% 
    right_join(mirx_labeled, by = c("read" = "sequence")) %>% 
    filter(is.na(mi_rna)) %>% .[["id"]]


# already missing
setdiff(mirx_labeled$id,clean_data$mi_rna) %>% .[grepl("hsa-",.[])]

# rank each sequence by each miRNA, and sample
parsed = 
    clean_data %>%
    filter(mi_rna == id | is.na(id)) %>% # only allow families where annotation is equal to expected 31312 rows
    annotate %>%
    group_by(mi_rna, sample) %>% # remove miRNAs that have more than one spike in in the family
    mutate(any_in_mirx = sum(grepl("hsa-", id))) %>% 
    ungroup() # 31312 seq,mir rows

# already missing
setdiff(mirx_labeled$id,parsed$mi_rna) %>% .[grepl("hsa-",.[])]

equimolar = parsed %>% filter(any_in_mirx == 1) %>% # 25665 (seq, mir rows), only allow families where the reference appears once
    left_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis

# We missed these families
setdiff(mirx_labeled$id,equimolar$mi_rna) %>% .[grepl("hsa-",.[])]

## mirge20 equimolar
gff_mirge = read_tsv("tools/mirge20/mirtop/expression_counts.tsv.gz") %>% 
    janitor::clean_names()
srr_to_name = read_csv("samples_bcbio_prepare.csv") %>% 
    mutate(samplename = tolower(samplename))

counts = gff_mirge[,13:ncol(gff_mirge)] %>% 
    as.matrix()
dds = DGEList(counts)
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = FALSE)
normalized = cbind(gff_mirge[,"read"],
                   counts) %>% 
    gather(sample, normalized,  -read)

# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
clean_data = gff_mirge %>% 
    left_join(mirx_labeled, by = c("read" = "sequence")) %>% # 31493 rows
    distinct(read, mi_rna, .keep_all = TRUE) # 31443 rows, need to identify what sequences are lost here

gff_mirge %>% 
    right_join(mirx_labeled, by = c("read" = "sequence")) %>% 
    filter(is.na(mi_rna)) %>% .[["id"]]

# already missing
setdiff(mirx_labeled$id,clean_data$mi_rna) %>% .[grepl("hsa-",.[])]

# rank each sequence by each miRNA, and sample
parsed = 
    clean_data %>%
    filter(mi_rna == id | is.na(id)) %>% # only allow families where annotation is equal to expected 31312 rows
    annotate %>%
    group_by(mi_rna, sample) %>% # remove miRNAs that have more than one spike in in the family
    mutate(any_in_mirx = sum(grepl("hsa-", id))) %>% 
    ungroup() # 31312 seq,mir rows

# already missing
setdiff(mirx_labeled$id,parsed$mi_rna) %>% .[grepl("hsa-",.[])]

equimolar_mirge = parsed %>% filter(any_in_mirx == 1) %>% # 25665 (seq, mir rows), only allow families where the reference appears once
    left_join(normalized, by = c("sample", "read")) %>% 
    mutate(sample = gsub("_cut.*$", "", sample)) %>% 
    left_join(srr_to_name, by = c("sample" = "samplename")) %>% 
    mutate(sample = description) %>% 
    select(-description) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis

# We missed these families
setdiff(mirx_labeled$id,equimolar_mirge$mi_rna) %>% .[grepl("hsa-",.[])]


# plasma
gffp = read_tsv("../tewari/tools/bcbio/mirtop/expression_counts.tsv.gz") %>% 
    janitor::clean_names()
meta_pilot = read_csv("../tewari/meta_pilot.csv")
id = meta_pilot$fixed_name
names(id) = meta_pilot$fixed_name
meta_pilot$sample_clean = names(janitor::clean_names(id))

counts = gffp[,13:ncol(gffp)] %>% 
    as.matrix()
dds = DGEList(counts)
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = FALSE)
normalized = cbind(gffp[,"read"],
                   counts) %>% 
    gather(sample, normalized,  -read) %>% 
    mutate(sample = gsub("_mirbase_ready", "", sample))
    
# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
plasma = gffp %>% 
    distinct(read, mi_rna, .keep_all = TRUE) %>%
    mutate(id = mi_rna) %>% 
    annotate %>%
    mutate(sample = gsub("_mirbase_ready", "", sample)) %>%
    inner_join(meta_pilot, by = c("sample" = "sample_clean")) %>%
    left_join(normalized, by = c("sample", "read")) %>% 
    rename(protocol=lib_method_simple) %>% 
    unite("short", c("protocol", "lab", "replicate"), remove = FALSE) %>% 
    distinct() %>% 
    analysis



# custom
gffc = read_tsv("custom_sequences/expression_counts.tsv.gz") %>% 
    janitor::clean_names() %>%
    select(-x4n_nex_tflex_lab7_synth_eq_clean) %>% 
    filter(nchar(iso_snp_nt) < 5)
counts = gffc[,13:ncol(gffc)] %>% 
    as.matrix()
dds = DGEList(counts)
dds = calcNormFactors(dds)
counts = cpm(dds, normalized.lib.sizes = FALSE)
normalized = cbind(gffc[,"read"],
                   counts) %>% 
    gather(sample, normalized,  -read)

custom = gffc %>% 
    distinct(read, mi_rna, .keep_all = TRUE) %>%
    mutate(id = mi_rna) %>% 
    annotate %>% 
    left_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis

save(custom, plasma, equimolar, equimolar_mirge, file = "data/data_gff.rda")

