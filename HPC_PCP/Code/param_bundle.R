library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

sim_out <- tibble()
xx <- 1

for (i in 1:42000) {
  if (file.exists(paste0("~/pcp_param/param_out_", i, ".RDA"))) {
  load(paste0("~/pcp_param/param_out_", i, ".RDA"))
  sim_out[xx,1:12] <- all_out
  xx <- xx + 1
  }
}

sim_out <- sim_out %>% 
  mutate_at(vars(2:3), as.factor) %>% 
  dplyr::select(-X, -L, -S, -L_out, -S_out)

save(sim_out, file = "param_out_all.RDA")
