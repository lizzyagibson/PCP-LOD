library(haven)
library(tidyverse)
library(janitor)
library(pcpr)
library(here)
library(PCPhelpers)

raw = read_xpt("./Data/L28POC_B.XPT")
names(raw)

all_nhanes = raw %>% 
      dplyr::select(-(2:3), -(148:165), -(169:195)) %>% 
      dplyr::select(-grep("[0-9]$", names(.)))

all_nhanes %>% drop_na
seqn = all_nhanes$SEQN

# Variables that end in "la" are the POP concentrations
conc_nhanes = all_nhanes %>% dplyr::select(grep("LA$", colnames(.)))
# Variables that end in "lc" are the <LOD labels
orig_label_nhanes = all_nhanes %>% dplyr::select(grep("LC$", colnames(.)))

# 1 == <LOD, change to -1
# 0 == >LOD, change to 1
replace_label = orig_label_nhanes
replace_label[replace_label == 1] = -1
replace_label[replace_label == 0] = 1

# element-wise multiply to set conc <LOD to -1
nhanes = (conc_nhanes * replace_label) %>% as_tibble()
nhanes[nhanes < 0] = -1

nhanes = cbind(seqn, nhanes)
label_nhanes = cbind(seqn, orig_label_nhanes)

pop_label_groups <- function(names) {
  names = case_when(grepl("(^D|d$|TCDD)", names) ~ "Dioxins", 
                    grepl("(^F|f$)", names) ~ "Furans",
                    grepl("(126|169|hxc|pcb|HXC|PCBPCB)", names) ~ "Non-Ortho PCBs",
                    grepl("(105|118|156|157|167|198)", names) ~ "Mono-Ortho PCBs",
                    TRUE ~ "Non-Dioxin-like PCBs")
  names = as_factor(names)
  names
}

pop_rename = function(name) {
  name = str_to_upper(name)
  name = case_when(grepl('074', name) ~ 'PCB 74',
                   grepl('099', name) ~ 'PCB 99',
                   grepl('118', name) ~ 'PCB 118',
                   grepl('138', name) ~ 'PCB 138',
                   grepl('153', name) ~ 'PCB 153',
                   grepl('170', name) ~ 'PCB 170',
                   grepl('180', name) ~ 'PCB 180',
                   grepl('187', name) ~ 'PCB 187',
                   grepl('194', name) ~ 'PCB 194',
                   grepl('D03', name) ~   '1,2,3,6,7,8-hxcdd',
                   grepl('D05', name) ~ '1,2,3,4,6,7,8-hpcdd',
                   grepl('D07', name) ~'1,2,3,4,6,7,8,9-ocdd',
                   grepl('F03', name) ~     '2,3,4,7,8-pncdf',
                   grepl('F04', name) ~   '1,2,3,4,7,8-hxcdf',
                   grepl('F05', name) ~   '1,2,3,6,7,8-hxcdf',
                   grepl('F08', name) ~ '1,2,3,4,6,7,8-hxcdf',
                   grepl('(hxc|HXC)', name) ~ 'PCB 169',
                   grepl('(LBDPCBLC|PCBPCB|LBXPCBLA)', name) ~ 'PCB 126', 
                   grepl("105", name) ~ "PCB 105",    
                   grepl("156", name) ~ "PCB 156", 
                   grepl("157", name) ~ "PCB 157", 
                   grepl("167", name) ~ "PCB 167",
                   grepl("D01", name) ~ "Dioxin 1",  
                   grepl("D04", name) ~ "Dioxin 4",   
                   grepl("F01", name) ~  "Furan 1",
                   grepl("F02", name) ~  "Furan 2",  
                   grepl("F04", name) ~  "Furan 4",  
                   grepl("F06", name) ~  "Furan 6", 
                   grepl("F07", name) ~  "Furan 7",   
                   grepl("TCD", name) ~     "TCDD",
                   grepl("189", name) ~  "PCB 189",  
                   grepl("196", name) ~  "PCB 196", 
                   grepl("199", name) ~  "PCB 199", 
                   grepl("D02", name) ~ "Dioxin 2", 
                   grepl("F09", name) ~  "Furan 9"  )
  name = as_factor(name)
  name
}

col_rename = function(pops) {
  colnames(pops) <- str_to_upper(colnames(pops))
  colnames(pops) <- str_sub(colnames(pops), 1, 6)
  colnames(pops) <- str_replace(colnames(pops), "(LBXD|LBDD)", "D")
  colnames(pops) <- str_replace(colnames(pops), "(LBXF|LBDF)", "F")
  colnames(pops) <- str_replace(colnames(pops), "TCD", "TCDD")
  colnames(pops) <- str_replace(colnames(pops), "(LBX|LBD)", "PCB")
  if ("PCBPCB" %in% colnames(pops)) {
    colnames(pops) <- str_replace(colnames(pops), "HXC", "169")
    colnames(pops) <- str_replace(colnames(pops), "PCBPCB", "PCB126")
  }
  return(pops)
}

# this takes proportion detected
# Some values were NA, no result, but not <LOD
# so go back to original labels
# 1 == <LOD
# 0 == >LOD
# ignore NA
prop <- function (x) {1-(sum(x, na.rm = TRUE)/length(x[!is.na(x)]))}

label_nhanes %>% 
  summarize_all(prop) %>% 
  pivot_longer(LBD052LC:LBDPCBLC) %>%
  mutate(name = str_to_upper(name),
         #name = pop_rename(name),
         name = fct_reorder(name, value)) %>% 
  ggplot(aes(x = name, y = value)) +
  geom_segment( aes(x=name, xend=name, y=0, yend=value), color="grey") +
  geom_point( color="orange", size=2) +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "% > LOD", n = "")

detected = label_nhanes %>%
  summarize_all(prop) %>% 
  pivot_longer(lbd074lc:lbdf09lc) %>%
  mutate(name = str_to_upper(name),
         name = pop_rename(name),
         name = fct_reorder(name, value)) %>% rename(Chemicals = name, Detected = value) %>% 
  arrange(desc(Detected))

# x is a proportion
# take all variables > x and preprocess them
# scale them and save their lods
process_pops = function(x) {
  # Select pops with detection above a certain limit
  props <- label_nhanes %>% dplyr::select(-seqn) %>%  summarize_all(prop)
  lod_names  = props %>% select_if(~. > x) %>% names()
  conc_names = lod_names %>% str_sub(., 4, 6) %>% str_c("lbx", ., "la")  
  
  pops <- nhanes %>% dplyr::select(seqn, all_of(conc_names)) %>% na.omit(.)
  seqn_ids = pops$seqn
  pops = col_rename(pops) %>% dplyr::select(-SEQN)
  
  # Make matrix of 0/1
  lods <- label_nhanes %>% 
    select(all_of(lod_names)) %>% 
    na.omit()
  
  # Matrix of all values (with CDC imputed values)
  imputed <- all_nhanes %>% 
    select(all_of(conc_names)) %>% 
    na.omit(.) %>% 
    col_rename(.) %>% 
    as.matrix()
  
  # Element-wise multiplication by 1/0
  # Keep values <LOD
  # Values >LOD == 0
  # Multiply by sqrt(2) to get back LOD
  lod_matrix <- (lods * imputed) * sqrt(2)
  lod_matrix <- as.matrix(lod_matrix)
  lod_matrix <- col_rename(lod_matrix)
  #summary(lod_matrix)
  
  # Need to scale POPs because they have super different ranges. 
  # Need to get rid of values less than LOD to scale and then add them back.
  
  # make <LOD NA so they dont affect the scaling
  pops[pops < 0] <- NA
  
  # Get stand dev of values > LDO
  denoms = apply(pops, 2, function(a) sd(a, na.rm = T))
  pops_scaled = apply(pops, 2, function(a) a/sd(a, na.rm = T))
  
  # make <LOD negative again
  pops_scaled[is.na(pops_scaled)] <- -1
  pops_scaled = as.matrix(pops_scaled)
  
  # Also want to scale the LODS!
  # Scale lod matrix by stand dev of measurements
  lod_matrix <- lod_matrix/denoms
  
  # and scale imputed
  imputed_scaled = apply(imputed, 2, function(a) a/sd(a, na.rm = T))
  # everything is scaled
  return(list(scaled_data = pops_scaled, lods = lod_matrix, sqrt2_data = imputed_scaled, seqn_ids = seqn_ids))
}
