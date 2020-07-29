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
