library(R.matlab)
library(tidyverse)
source("./PCP/stable_pcp_alternating.R")

# Load data
mixture <- readMat("./mixtures_data.mat")

mixture_data <- as.data.frame(mixture) %>% as.tibble() %>% select(Al, As, Ba, bc, Br, Ca, Cl,
                                                                  Cr, Cu, Fe, K,  Mn,  Ni,  Pb,  S,  Se,  Si,
                                                                  Ti,  V, Zn)

mix <- mixture_data[complete.cases(mixture_data),]
X <- scale(mix, center = TRUE, scale = TRUE)

m <- nrow(X)
n <- ncol(X)

lambda = 1/sqrt(m)

mixture_out <- stable_pcp_alternating(X, 4/sqrt(m), 10)

mixture_S <- mixture_out$S
mixture_L <- mixture_out$L

# Plot
mixture_S %>% as.tibble() %>% 
  mutate(id = 1:nrow(mixture_S)) %>% select(id, everything()) %>% 
  gather(key = exposure, value = value, Al:Zn) %>%
  ggplot(aes(x = exposure, y = id)) +
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "navy", mid = "navy", high = "yellow", 
                       na.value = "transparent") +
  labs(x = "Exposure", y = "Date", title = "Sparse matrix of rare events", legend = "Magnitude") + 
  theme_classic()

#Plot 1a (specific dates)
dateLabels = c('5/28/10','5/29/10','5/30/10','5/31/10','6/1/10','6/2/10','6/3/10','6/4/10','6/5/10')

figure1 <- cbind(dateLabels, mixture_S[(2417-3):(2417+5),]) %>% as.tibble() %>% 
  gather(key = exposure, value = value, Al:Zn) %>%
  mutate(value = as.numeric(value)) %>% 
  ggplot(aes(x = exposure, y = dateLabels)) +
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "navy", mid = "navy", high = "yellow",
                       na.value = "transparent") +
  labs(x = "Exposure", y = "Date", title = "Sparse matrix of rare events", legend = "Magnitude") + 
  theme_classic()

png("figure1.png", width = 2000, height = 1200, res = 300)
figure1
dev.off()

# Plot 2
figure2 <- as.tibble(svd(mixture_L)$v) %>% 
  mutate(id = c("Al", "As", "Ba", "bc", "Br", "Ca", "Cl",
                "Cr", "Cu", "Fe", "K",  "Mn",  "Ni",  "Pb",  "S",  "Se",  "Si",
                "Ti",  "V", "Zn")) %>% 
  select(id, everything()) %>% 
  gather(key = singular_vector, value = magnitude, V1:V20) %>%
  filter(singular_vector %in% c("V1", "V2", "V3", "V4", "V5")) %>% 
  ggplot(aes(x = id, y = magnitude)) + geom_point(color = "blue") + 
  geom_segment(aes(xend = id, yend = 0), color = "blue") +
  facet_grid(. ~ singular_vector) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "red") +
  theme_bw() + labs(x = "", y = "Magnitude", title = "First 5 right singular vectors of svd(L)")

png("figure2.png", width = 5000, height = 1500, res = 300)
figure2
dev.off()


