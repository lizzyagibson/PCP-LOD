% demo_mixtures
clear;

addpath('/Users/lizzy/OneDrive - cumc.columbia.edu/Principal.Component.Pursuit');
load('/Users/lizzy/OneDrive - cumc.columbia.edu/Principal.Component.Pursuit/Data/mixtures_data.mat');

%X = [pm25 pm1 Al As Ba bc Br Ca Cl Cr Cu Fe K Mn Ni Pb S Se Si Ti V Zn];
X = [Al As Ba bc Br Ca Cl Cr Cu Fe K  Mn  Ni  Pb  S  Se  Si Ti  V Zn];
n = [ 1  2  3  4  5  6  7  8  9 10 11 12  13  14  15 16  17 18 19 20];
% [] makes a vector, take columns from data set, put into matrix as columns

numMissingPerRow = sum( isnan(X), 2 ); 
%get rid of rows with NANs
goodRows = find( numMissingPerRow == 0 ); 
% good rows without missing data

X = X(goodRows,:); 
%semicolon means it doesnt output the results

[m,n] = size(X);
% m and n become the number of rows and columns

Xmissing = X;
Xmissing(1:1000,1:5) = NaN;

lambda = 1/sqrt(m); 
mu = 1;

%% Run models
% [L,S] = root_pcp_rank_r_with_missing(Xmissing, lambda, 1, 5); 
% norm(L, "Fro")
% norm(S, "Fro")

[L,S] = root_pcp_rank_r_nonnegL_with_missing(Xmissing, lambda, 1, 5); 
norm(L, "Fro")
norm(S, "Fro")


% [Llod, Slod, losslod] = pcp_lod(X, lambda, mu, 0); 
% %Converged in 28 iterations.
% norm(Llod, "Fro")
% norm(Slod, "Fro")

% [Lroot,Sroot] = root_pcp_rank_r_nonnegL(X, lambda, 1, 5); 
% %Converged in 837 iterations.
% norm(Lroot, "Fro")
% norm(Sroot, "Fro")

% [L,S] = root_pcp_nonnegL(X, lambda, mu);
% %Converged in 646 iterations.
% norm(L, "Fro")
% norm(S, "Fro")

% [Lnan,Snan] = root_pcp_with_nan(Xmissing, lambda, mu);
% % %Converged in 915 iterations.
% norm(Lnan, "Fro")
% norm(Snan, "Fro")

% [Lnan_nn,Snan_nn] = root_pcp_with_nan_nonnegL(Xmissing, lambda, mu);
% %Converged in 646 iterations.
% norm(Lnan_nn, "Fro")
% norm(Snan_nn, "Fro")

% [Lnan_nn_lod,Snan_nn_lod] = root_pcp_with_nan_nonnegL_LOD(Xmissing, lambda, mu, 0);
% %Converged in 632 iterations.
% norm(Lnan_nn_lod, "Fro")
% norm(Snan_nn_lod, "Fro")

