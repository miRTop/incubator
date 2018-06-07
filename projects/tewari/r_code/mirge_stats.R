setwd(here::here())
library(tidyverse)
library(ggplot2)
theme_set(
    theme_light(base_size = 9))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

meta_pilot = read_csv("meta_pilot.csv")

## needs to get name inside the file and name of the file to map to original
keys = read_csv("tools/mirge/sample_fn_key.txt", col_names = c("fn", "inside")) %>% 
    mutate(fn = gsub("_isomiRs.gff", "", fn))
mirge_stats = read_csv("tools/mirge/stats/mirtop_stats.txt", skip = 1) %>% 
    inner_join(keys, by = c("sample" = "inside")) %>% 
    inner_join(meta_pilot,
               by = c("fn" = "fixed_name")) %>% 
    mutate(sample = fn)

save_stats(mirge_stats, "miRge")