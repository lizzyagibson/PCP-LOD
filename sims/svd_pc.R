load("./sims/pcp_out_all.rda")

# SVD ####
svd_out = pcp_out %>% 
  mutate(svd_chem = map(chem, svd),
         svd_pca  = map(lod_sqrt2_mat, svd),
         L        = map(pcp_out, function(x) x$L),
         svd_pcp  = map(L, svd),
         svd_chem_left  = map(svd_chem, function(x) x$u),
         svd_chem_right = map(svd_chem, function(x) x$v),
         svd_pca_left  = map(svd_pca, function(x) x$u),
         svd_pca_right = map(svd_pca, function(x) x$v),
         svd_pcp_left  = map(svd_pcp, function(x) x$u),
         svd_pcp_right = map(svd_pcp, function(x) x$v))

class(svd_out$svd_pcp_left[[1]])

svd_re1 = svd_out %>% 
  mutate(svd_pca_left  = map2(svd_chem_left, svd_pca_left, factor_correspondence))
svd_re01 = svd_re1 %>% dplyr::select(1:18, svd_pca_left)

svd_re2 = svd_out %>% 
  mutate(svd_pcp_left  = map2(svd_chem_left, svd_pcp_left, factor_correspondence))
svd_re02 = svd_re2 %>% dplyr::select(1:18, svd_pcp_left)

svd_re3 = svd_out %>% 
  mutate(svd_pcp_right = map2(svd_chem_right, svd_pcp_right, factor_correspondence))
svd_re03 = svd_re3 %>% dplyr::select(1:18, svd_pcp_right)

svd_re4a = svd_out %>% 
  slice(1:100) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4b = svd_out %>% 
  slice(101:200) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4c = svd_out %>% 
  slice(201:300) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4d = svd_out %>% 
  slice(301:400) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4g = svd_out %>% 
  slice(401:424) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4f = svd_out %>% 
  slice(451:500) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4h = svd_out %>% 
  slice(425:450) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re4e = svd_out %>% 
  slice(501:600) %>% 
  mutate(svd_pca_right = map2(svd_chem_right, svd_pca_right, factor_correspondence))

svd_re = full_join(svd_re01, svd_re02) %>% full_join(., svd_re03)

svd_re04 = bind_rows(svd_re4a, svd_re4b, svd_re4c, svd_re4d, svd_re4g, svd_re4f, svd_re4h, svd_re4e) %>% 
  dplyr::select(1:18, svd_pca_right)

svd_re = full_join(svd_re, svd_re04)
# save(svd_re, file = "./sims/svd_out_re.rda")  

svd_metrics = svd_re %>% 
  mutate(svd_err_pca_left   = map2(svd_chem_left, svd_pca_left, function(x,y) norm(x-y,"F")/norm(x,"F")),
         svd_err_pcp_left   = map2(svd_chem_left, svd_pcp_left, function(x,y) norm(x-y,"F")/norm(x,"F")),
         svd_err_pca_right   = map2(svd_chem_right, svd_pca_right, function(x,y) norm(x-y,"F")/norm(x,"F")),
         svd_err_pcp_right   = map2(svd_chem_right, svd_pcp_right, function(x,y) norm(x-y,"F")/norm(x,"F"))) %>% 
  unnest(c(svd_err_pca_left, svd_err_pcp_left, svd_err_pca_right, svd_err_pcp_right)) %>% 
  group_by(chemicals, lim) %>% 
  summarize(pca_left = quantile(svd_err_pca_left), prop = seq(0,1,.25),
            pcp_left = quantile(svd_err_pcp_left),
            pca_right = quantile(svd_err_pca_right),
            pcp_right = quantile(svd_err_pcp_right))

svd_metrics %>% 
  pivot_longer(c(pca_left, pcp_left:pcp_right)) %>% 
  pivot_wider(names_from = prop,
              values_from = value) %>% print(n=24)
