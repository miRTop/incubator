library(readr)

fns = list.files("tables/dispersion", full.names = T)
fns = fns[!grepl("lib", fns)]

df = map(fns, read_csv) %>% bind_rows() %>% 
    mutate(pct_cat = factor(pct_cat, levels=c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")))

library(ggplot2)

lookat = c("add3p", "shift5p", "shift3p", "snp", "reference")

# general dispersion by lab-protocol
separate(df, group, into = c("lab", "protocol"), sep = " ") %>% 
    filter(n > 20) %>% 
    #unite(protocol, tool, col = "gg_group", remove = F) %>% 
    ggplot(aes(lab,disp,color=protocol)) +
    geom_boxplot() +
    ggsave("figures/dispersion/general_lab_protocol.pdf")
    
# dispersion by lab - one tool - by isomir
separate(df, group, into = c("lab", "protocol"), sep = " ") %>% 
    filter(tool == "bcbio", n > 20, iso  %in% lookat, protocol=="NEBNext") %>% 
    #unite(protocol, tool, col = "gg_group", remove = F) %>% 
    ggplot(aes(exp<5, disp, fill = pct_cat)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(lab~iso) +
    theme(axis.text = element_text(size = 9),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    ggtitle("bcbio - NEBNext")
separate(df, group, into = c("lab", "protocol"), sep = " ") %>% 
    filter(tool == "bcbio", n > 20, iso  %in% lookat, protocol=="TrueSeq") %>% 
    #unite(protocol, tool, col = "gg_group", remove = F) %>% 
    ggplot(aes(exp<5, disp, fill = pct_cat)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(lab~iso) +
    theme(axis.text = element_text(size = 9),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    ggtitle("bcbio - TrueSeq")

# dispersion by isomir type and importance
separate(df, group, into = c("lab", "protocol"), sep = " ") %>% 
    filter( n > 20, iso  %in% lookat) %>% 
    unite(protocol, tool, col = "gg_group", remove = F) %>% 
    ggplot(aes(x = pct_cat, y = disp, color = tool, linetype = protocol,
               group = gg_group)) +
    geom_line() +
    facet_grid(lab~iso) +
    theme(axis.text = element_text(size = 9),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    ggsave("figures/dispersion/specific_isomir_lab_tool_protocol.pdf")

separate(df, group, into = c("lab", "protocol"), sep = " ") %>% 
    filter( n > 20, iso  %in% lookat) %>% 
    unite(protocol, tool, col = "gg_group", remove = F) %>% 
    ggplot(aes(x = pct_cat, y = log10(n), color = tool, linetype = protocol,
               group = gg_group)) +
    geom_line() +
    facet_grid(lab~iso) +
    theme(axis.text = element_text(size = 9),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))
separate(df, group, into = c("lab", "protocol"), sep = " ") %>% 
    filter( n > 20, iso  %in% lookat) %>% 
    unite(protocol, tool, col = "gg_group", remove = F) %>% 
ggplot(aes(x = pct_cat, y = (exp), color = tool, linetype = protocol,
           group = gg_group)) +
    geom_line() +
    facet_grid(lab~iso) +
    theme(axis.text = element_text(size = 9),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))

# no correlation between number and dispersion
ggplot(df, aes(x = log10(n), y = disp, fill = pct_cat)) +
    geom_point() +
    facet_grid(tool~iso) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))

# no correlation between exp and dispersion
ggplot(df, aes(x = exp, y = disp, fill = pct_cat)) +
    geom_point() +
    facet_grid(tool~iso) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))
