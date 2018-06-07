setwd(here::here())
library(tidyverse)
library(ggplot2)
theme_set(
    theme_light(base_size = 9))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

meta_pilot = read_csv("meta_pilot.csv")

stats = read_csv("tools/srnabench/mirtop/mirtop_stats.txt", comment = "#") %>% 
    inner_join(meta_pilot,
               by = c("sample" = "fixed_name"))

save_stats(stats, "sRNAbench")
