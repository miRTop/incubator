setwd(here::here())
library(tidyverse)
library(ggplot2)
library(cowplot)
source("r_code/functions.R")
theme_set(theme_bw(base_size = 9))

meta_pilot = read_csv("meta_pilot.csv")

ggplot(meta_pilot, aes(lab, fill = lib_method_simple)) +
    geom_bar(position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    ggsave("figures/pilot.png")
