library(tidyverse)
library(seqinr)
setwd(here::here()) # set working directory in project folder
source("r_code/seqThreshold.R")
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
    group_by(read) %>% 
    mutate(hits=n()) %>% filter(hits==1) %>% 
    mutate(fix_add = iso_add ==  -1 * iso_3p & iso_add != 0,
           iso_add = ifelse(fix_add, 0, iso_add),
           iso_3p = ifelse(fix_add, 0, iso_3p),
           iso_snp = ifelse(fix_add, 1, iso_snp)) %>% 
    mutate(is_snp = ifelse(iso_snp != 0, "snp", ""),
           is_add = ifelse(iso_add != 0, "add3p", ""),
           is_5p = ifelse(iso_5p != 0, "shift5p", ""),
           is_3p = ifelse(iso_3p != 0, "shift3p", ""),
           ) %>% 
    mutate(is_loss_5p = ifelse(iso_5p < 0, "loss_5p", ""),
           is_loss_5p = ifelse(iso_5p > 0, "gain_5p", is_loss_5p),
           is_loss_3p = ifelse(iso_3p < 0, "loss_3p", ""),
           is_loss_3p = ifelse(iso_3p > 0, "gain_3p", is_loss_3p)) %>%
    unite("iso_shift_nt", iso_5p, iso_3p, sep = " ", remove = FALSE) %>%
    unite("iso_nt", iso_snp, iso_add, iso_5p, iso_3p, sep = "_") %>%
    unite("iso_detail", iso_add_nt, iso_5p_nt, iso_3p_nt, sep = "_") %>%
    unite("iso", is_snp, is_add, is_5p, is_3p, sep = " ") %>%
    mutate(iso = trimws(iso),
           iso = gsub("[ ]+", " ", iso),
           iso = ifelse(iso=="", "reference",iso)) %>% 
    unite("iso_shift", is_loss_5p, is_loss_3p, sep =  " ") %>%
    mutate(iso_shift = trimws(iso_shift),
           iso_shift = ifelse(iso_shift=="", "no-loss",iso_shift)) %>% 
    select(-uid, -variant, -iso_snp_nt, -fix_add, -hits ) %>% 
    ungroup %>% 
    gather(sample, value, -mi_rna, -iso, -iso_nt, -iso_shift, -iso_shift_nt, -id, -read, -iso_detail) %>% 
    filter(value>0)

# quantify the times sequences are detected and % of importance
analysis = . %>% 
    group_by(sample, mi_rna) %>% 
    filter(value > 0) %>% 
    arrange(sample, mi_rna, desc(value)) %>%
    mutate(rank = 1:n(),
           total = sum(value),
           pct = value / total * 100) %>% 
    ungroup() %>% 
    group_by(read) %>% 
    mutate(reproducible_protocol = paste(unique(protocol), collapse = " "),
           reproducible_lab = length(unique(lab))) %>%
    group_by(mi_rna, sample) %>%
    mutate(ref_is_1 = length(mi_rna[rank == 1 & iso == "reference"])) %>% 
    ungroup %>% 
    mutate(pct_cat = cut(pct,
                          breaks = c(-1, 0.1, 1, 5, 10, 20, 50, 101),
                          labels = c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")),
           mir_coverage = cut(total,
                              breaks = c(-1, 1, 10, 100, 1000, 1e30),
                              labels = c("<1", "1-10", "10-100", "100-1000", ">1000"))) %>% 
    ungroup()

norm = function(fn, cache=NULL){
    counts = fn[,13:ncol(fn)] %>% 
        as.matrix()
    counts = counts[,colSums(counts) > 10]
    if (!is.null(cache)){
        if (file.exists(cache)){
            load(cache)
        }else{
            t = apply(counts, 2, function(c){
                thresholdSeq(c[c>0])[["threshold"]]
            })
            threshold = data.frame(sample = names(t), 
                                   threshold = t, stringsAsFactors = F)
        }

    }
    save(threshold, file = cache)
    dds = DGEList(counts)
    dds = calcNormFactors(dds)
    counts = cpm(dds, normalized.lib.sizes = FALSE)
    normalized = cbind(fn[,"read"],
                       counts) %>% 
        gather(sample, normalized,  -read) %>% 
        left_join(threshold, by = "sample")
    
}

# equimolar
gff = read_tsv("tools/bcbio/mirtop/expression_counts.tsv.gz") %>% 
    janitor::clean_names()

normalized = norm(gff, "data/th_gff.rda")
# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
equimolar = gff %>% 
    left_join(mirx_labeled, by = c("read" = "sequence"))  %>% # mirx_labeled defined in line 14, it contains the mirxplor sequences and the updates names
    filter(mi_rna == id | is.na(id)) %>% # only allow families where annotation is equal to expected or a variant if is NA
    annotate %>%
    group_by(mi_rna, sample) %>% # remove miRNAs that have more than one spike in in the family
    mutate(any_in_mirx = sum(grepl("hsa-", id)),
           mirx_counts = length(unique(id[!is.na(id)]))) %>% 
    ungroup() %>% 
    filter(any_in_mirx == 1, mirx_counts == 1) %>%  # where the miRNA matches only once with the mirxplor sequences
    left_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis

# We missed these families
# setdiff(mirx_labeled$id,equimolar$mi_rna) %>% .[grepl("hsa-",.[])]

## mirge20 equimolar
gff_mirge = read_tsv("tools/mirge20/mirtop/expression_counts.tsv.gz") %>% 
    janitor::clean_names()
srr_to_name = read_csv("config/samples_bcbio_prepare.csv") %>% 
    mutate(samplename = tolower(samplename))

normalized = norm(gff_mirge, "data/th_gffmirge.rda")

equimolar_mirge = gff_mirge %>% 
    left_join(mirx_labeled, by = c("read" = "sequence")) %>% 
    filter(mi_rna == id | is.na(id)) %>% # only allow families where annotation is equal to expected in mirxplor (same id in the analysis and mirxplor file)
    annotate %>%
    group_by(mi_rna, sample) %>% # remove miRNAs that have more than one spike in in the family
    mutate(any_in_mirx = sum(grepl("hsa-", id)), # numbers of human mirnas matching the mirxplor sample. Ideally one per miRNA if the mirna has not similarity with others
       mirx_counts = length(unique(id[!is.na(id)]))) %>% # numbers of mirnas matching the mirxplor sample. Ideally one per miRNA...
    ungroup() %>% 
    filter(any_in_mirx == 1, mirx_counts == 1) %>% 
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

# equimolar razer3
gffr = read_tsv("tools/razer3/expression_counts.tsv.gz") %>% 
    janitor::clean_names() %>% 
    filter(nchar(iso_snp_nt) < 5)

normalized = norm(gffr, "data/th_gffr.rda")

# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
equimolar_razer3 = gffr %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
    left_join(mirx_labeled, by = c("read" = "sequence")) %>% 
    annotate %>% # fix naming
    group_by(mi_rna, sample) %>% # remove miRNAs that have more than one spike in in the family or not the correct one
    mutate(any_in_mirx = sum(grepl("hsa-", id)),
           mirx_counts = length(unique(id[!is.na(id)]))) %>% 
    ungroup() %>% 
    filter(any_in_mirx == 1, mirx_counts == 1) %>% 
    inner_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"), # fix metadata
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis # categorize isomirs


# dsrg
gffdsrg = read_tsv("tools/dsrg/expression_counts.tsv.gz") %>% 
    janitor::clean_names()  %>% 
    set_names(gsub("_mirbase_ready", "", names(.)))

normalized = norm(gffdsrg, "data/th_gffdsrg.rda")

# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
dsrg = gffdsrg %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
    set_names(gsub("_mur_d", "_murd", names(.))) %>%
    .[,c(1:12, grep("_mur_", names(.)))] %>% 
    .[rowSums(as.matrix(.[,13:ncol(.)]))>0,] %>% 
    left_join(mirx_labeled, by = c("read" = "sequence")) %>% 
    filter(mi_rna == id | is.na(id)) %>% 
    annotate %>%
    group_by(mi_rna, sample) %>% 
    mutate(any_in_mirx = sum(grepl("hsa-", id)),
           mirx_counts = length(unique(id[!is.na(id)]))) %>% 
    ungroup() %>% 
    filter(any_in_mirx == 1, mirx_counts == 1) %>% 
    mutate(sample = gsub("mur_d", "murd", sample)) %>% 
    inner_join(normalized, by = c("sample", "read")) %>% 
    separate(sample, remove = F, into = c("protocol", "source", "lab", "index", "snumber")) %>% 
    unite("short", c("protocol", "source", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis


# plasma
gffp = read_tsv("../tewari/tools/bcbio/mirtop/expression_counts.tsv.gz") %>% 
    janitor::clean_names()
meta_pilot = read_csv("../tewari/meta_pilot.csv")
id = meta_pilot$fixed_name
names(id) = meta_pilot$fixed_name
meta_pilot$sample_clean = names(janitor::clean_names(id))

normalized = norm(gffp, "data/th_gffp.rda") %>% 
    mutate(sample = gsub("_mirbase_ready", "", sample))
    
# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
plasma = gffp %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
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
gffc = read_tsv("tools/custom_sequences/expression_counts.tsv.gz") %>% 
    janitor::clean_names() %>%
    select(-x4n_nex_tflex_lab7_synth_eq_clean) %>% 
    filter(nchar(iso_snp_nt) < 5)
normalized = norm(gffc, "data/th_gffc.rda")


custom = gffc %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
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


# vandijk
gffv = read_tsv("tools/vandijk/expression_counts.tsv.gz") %>% 
    janitor::clean_names() %>%
    filter(nchar(iso_snp_nt) < 5)
normalized = norm(gffv, "data/th_gffv.rda")


vandijk = gffv %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
    mutate(id = mi_rna) %>% 
    annotate %>% 
    inner_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"rep[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis %>% 
    mutate(replicate=stringr::str_extract(short, "[0-9]{1,2$"))


# narrykim
gffnk = read_tsv("tools/narrykim/spikeins/mirtop.tsv.gz") %>% 
    janitor::clean_names() %>%
    filter(nchar(iso_snp_nt) < 5)
normalized = norm(gffnk, "data/th_gffnk.rda")


narrykim = gffnk %>% 
    rename(iso_add=iso_add3p,
           iso_add_nt=iso_add3p_nt) %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
    mutate(id = mi_rna) %>% 
    annotate %>% 
    inner_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"^\\w{1,4}_"),
           protocol=stringr::str_extract(sample, "_\\w{1,2}_"),
           index =stringr::str_extract(sample, "[0-9]$")) %>% 
    mutate(short=sample) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis %>% 
    mutate(replicate=stringr::str_extract(short, "[0-9]$"))


# narrykim2
gffnk = read_tsv("tools/narrykim/human/mirtop.tsv.gz") %>% 
    janitor::clean_names() %>%
    filter(nchar(iso_snp_nt) < 5)
normalized = norm(gffnk, "data/th_gffnk2.rda")


narrykim_human = gffnk %>% 
    dplyr::rename(iso_add=iso_add3p,
           iso_add_nt=iso_add3p_nt) %>% 
    # distinct(read, mi_rna, .keep_all = TRUE) %>%
    mutate(id = mi_rna) %>% 
    annotate %>% 
    inner_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"^\\w{1,4}_"),
           protocol=stringr::str_extract(sample, "_\\w{1,2}_"),
           index =stringr::str_extract(sample, "[0-9]$")) %>% 
    mutate(short=sample) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis %>% 
    mutate(replicate=stringr::str_extract(short, "[0-9]$"))


save(equimolar_razer3, gffnk, vandijk, custom, plasma, equimolar, equimolar_mirge, dsrg, narrykim, narrykim_human,
     file = "data/data_gff.rda")
