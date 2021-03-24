data50 = table2array(readtable("/Users/lizzy/Principal.Component.Pursuit/experiments/lizzy_experiments/nhanes/data/nhanes_neg1.csv"));
delta50 = table2array(readtable("/Users/lizzy/Principal.Component.Pursuit/experiments/lizzy_experiments/nhanes/data/nhanes_lods.csv"));

[n, p] = size(data50);

lambda = 0.02459549; % still doesn't converge in 500000 iterations.
% lambda = 0.022; this one does converge in 22043 iterations.
% lambda = 0.028; this one does converge in  1240 iterations.

mu = sqrt(p/2);
r = 4;

[L,S] = root_pcp_ncvx_nan_nonnegL_LOD(data50, lambda, mu, r, delta50);

rank(L, 1e-04)
sum(sum(S >= -1e-04 & S <= 1e-04))/(n*p)
norm(L, 'fro')
norm(S, 'fro')