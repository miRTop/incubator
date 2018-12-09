library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

prepare =  . %>% filter(ref_is_1 == 1, nchar(iso_detail)<7) %>%
    dplyr::count(protocol, sample, pct_cat, iso_shift_nt, iso, iso_detail) %>% 
    group_by(protocol) %>% 
    mutate(sample2 = paste0(protocol, "_", 1:length(unique(sample)))) %>% 
    group_by(sample2, protocol) %>% 
    mutate(pct_shift = n/sum(n)*100)


background = equimolar_razer3 %>% 
    filter(ref_is_1 == 1, iso == "reference") %>% 
    distinct(read) %>% 
    mutate(nt1 = stringr::str_sub(read, 1, 1),
           ntlast = stringr::str_sub(read, nchar(read), nchar(read)))

bind_rows(
count(background, nt1) %>% 
    mutate(total = sum(n),
           pct = n/total * 100,
           position = "first") %>% 
    rename(nt = nt1),
count(background, ntlast) %>% 
    mutate(total = sum(n),
           pct = n/total * 100,
           position = "last") %>% 
    rename(nt = ntlast)
) %>% 
    ggplot(aes(nt, pct, fill = position)) +
    geom_bar(stat = "identity", position = "dodge") + 
    ggsave("results/06_shift_nts/06_background_shift_nts.pdf")


equimolar_razer3 %>% 
    filter(iso_shift_nt!="0", (iso=="shift5p" | iso=="shift3p"),
           iso_shift_nt  %in% c("-1 0", "0 -1")) %>% 
    prepare() %>% 
    mutate(nt = gsub("[_0]", "", iso_detail)) %>% 
    ggplot(aes(x = protocol, y = pct_shift, color = nt, fill = nt)) +
    geom_boxplot() +
    scale_fill_discrete() +
    facet_grid(pct_cat~iso_shift_nt) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.text = element_text(size = 14, color = "black")) +
    ylab("% of isomiRs") +
    ggsave("results/06_shift_nts/06_shift_nts.pdf", height = 9)
