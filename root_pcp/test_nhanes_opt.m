data50 = table2array(readtable("/Users/lizzy/experiments/lizzy_experiments/nhanes/data/nhanes_neg1.csv"));
delta50 = table2array(readtable("/Users/lizzy/experiments/lizzy_experiments/nhanes/data/nhanes_lods.csv"));

[n, p] = size(data50);

lambda = 0.02459549; % still doesn't converge in 500000 iterations.
% lambda = 0.022; % this one does converge in 22043 iterations.
% lambda = 0.028; % this one does converge in  1240 iterations.

mu = sqrt(p/2); % 3.2404
% mu = 3 % Converged in 1423 iterations.
% mu = 6

% lambda = 0.0076
% mu = 1
% Converged in 11939 iterations.

r = 4;

% mu = 131.7465
% lambda = 1

% lam2 = lambda/(lambda+mu)
% mu2 = mu/(lambda+mu)

[Lno,Sno] = root_pcp_ncvx_nan_nonnegL_LOD(data50, lambda, mu, r, delta50);

rank(Lno, 1e-04)
rank(Lyes, 1e-04)

sum(sum(Sno >= -1e-04 & Sno <= 1e-04))/(n*p)
sum(sum(Syes >= -1e-04 & Syes <= 1e-04))/(n*p)

norm(Lno, 'fro')
norm(Lyes, 'fro')


norm(Sno, 'fro')
norm(Syes, 'fro')

