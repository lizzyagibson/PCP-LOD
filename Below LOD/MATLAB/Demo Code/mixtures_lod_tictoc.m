% demo_mixtures
clear;

load('../../../Data/mixtures_data.mat');

X = [Al As Ba bc Br Ca Cl Cr Cu Fe K  Mn  Ni  Pb  S  Se  Si Ti  V Zn];
n = [ 1  2  3  4  5  6  7  8  9 10 11 12  13  14  15 16  17 18 19 20];

numMissingPerRow = sum( isnan(X), 2 ); 
goodRows = find( numMissingPerRow == 0 ); 
X = X(goodRows,:); 
[m,n] = size(X);

colStd = std(X);
X_std = X ./ repmat(colStd, m, 1);

lambda = 1/sqrt(m); 
%weight parameter for pcp, 1 / sqrt(number of rows)
mu = sqrt(n/(2*log(m*n)));

% % Run model reg
% tic()
% [L, S, loss] = pcp_lod(X, lambda, mu, 0); 
% toc()
% Elapsed time is 2.768229 seconds.
% Converges in 747

% % Run model standardized
tic()
[Lstd, Sstd, lossstd] = pcp_lod(X_std, lambda, mu, 0); 
toc()
% Elapsed time is 5.318982 seconds.
% Converges in 1905
