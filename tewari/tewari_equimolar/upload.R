library(rdrop2)
token <- readRDS("~/.droptoken.rds")
d = drop_acc(dtoken = token)
dropdir = "mirtop/synthetic"

library(r2dropSmart)
sync(".", dropdir, 
     pattern = "rda",
     token = token,
     share = T,
     dry = T)
#