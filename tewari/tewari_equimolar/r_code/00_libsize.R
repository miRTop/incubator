library(tidyverse)
library(ggplot2)
library(pheatmap)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")


equimolar_razer3 %>% 
    group_by(sample, protocol) %>% 
    summarise(n=n(),size=sum(value)) %>% 
    mutate(ytext = size + size*0.1,
           text = paste0(round(n/1000, digits = 0), "K")) %>% 
    ggplot(aes(sample, size, fill = protocol)) +
    geom_bar(stat = "identity") +
    ylab("reads") +
    geom_text(aes(sample, ytext, label = text), size = 3) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggsave("results/00_library_size/00_libsize_diversity_bar.pdf", width = 9)

equimolar_razer3 %>% 
    group_by(sample, protocol) %>% 
    summarise(n=n(),size=sum(value)) %>% 
    mutate(ytext = size + size*0.1,
           text = paste0(round(n/1000, digits = 0), "K")) %>% 
    ggplot(aes(size, n, color = protocol)) +
    scale_x_log10() +
    geom_point(size=3) +
    ggsave("results/00_library_size/00_libsize_diversity_scatter.pdf")
