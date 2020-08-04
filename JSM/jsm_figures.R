library(GGally)

n = 1000
r = 4
p = 20

U = matrix(runif(n*r), nrow = n, ncol = r)
V = matrix(runif(r*p), nrow = r, ncol = p)

L = U %*% V
apply(L, 2, sd)
sd(L)
var(as.vector(L))

as_tibble(L) %>% 
  mutate(id = 1:nrow(L)) %>% 
  pivot_longer(V1:V20) %>% 
  #filter(name == "V1") %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~name)

as_tibble(L) %>% 
  mutate(id = 1:nrow(L)) %>% 
  pivot_longer(V1:V20) %>% 
  #filter(name == "V1") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50, fill = "#1976D2") + # HTML color code to match latex
  theme_minimal() +
  labs(x = "Simulated Variable", y = "Count")

Z = matrix(rnorm(n*p, mean = 0, sd = .5*sd(L)), nrow = n, ncol = p)
X = L + Z
D = pmax(X, 0)
apply(D, 2, sd)
apply(D, 2, max)
sd(D)

quantile(L, probs = 0.015)

sim_hist <- as_tibble(L) %>% 
  mutate(id = 1:nrow(L)) %>% 
  pivot_longer(V1:V20) %>% 
  #filter(name == "V1") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50, fill = "#1976D2") + # HTML color code to match latex
  theme_minimal(base_size = 20) +
  labs(x = "Simulated Low Rank Matrix", y = "Count")

#library(GGally)
sim_corr <- ggcorr(as_tibble(D), limits = FALSE,
       hjust = 0.85, size = 3, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 20)

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/jsm_hist.pdf")
sim_hist
dev.off()

pdf("/Users/lizzy/Summer 2020/JSM PCP-LOD/figures/jsm_corr.pdf")
sim_corr
dev.off()

#########################################

sigma = 0.6 * sd(L)

Z <- mvtnorm::rmvnorm(n, mean = rep(0, 20), sigma = diag(0.1, 20, 20))
X = L + Z

D = pmax(L, 0)

Delta = quantile(D, 0.1)

D_minus1 = (D>=Delta)*D + (D<Delta)*(-1)
#D_minus1[D_minus1 == -1] = NA
head(D_minus1)
summary(D_minus1)

#library(pcpr)
#install.packages("formattable")
#library(formattable)

Dout <- pcp_lod(D_minus1, lambda = 1/sqrt(nrow(D)), mu = 38, LOD = Delta)

round(D[2:4, 2:4], 3)

round(D_minus1[2:4, 2:4], 3)

round(Dout$L[2:4, 2:4], 3)


