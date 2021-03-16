% demo_mixtures
clear;

addpath('/Users/lizzy/Principal.Component.Pursuit');
load('/Users/lizzy/Principal.Component.Pursuit/Data/mixtures_data.mat');

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

[m,p] = size(X);
% m and n become the number of rows and columns

Xmissing = X;
Xmissing(1:10,1:5) = NaN;

Xlod = X;
Xlod(10:20,6:10) = -1;

Xmissinglod = X;
Xmissinglod(1:10,1:5) = NaN;
Xmissinglod(10:20,6:10) = -1;

lambda = 1/sqrt(m); 
mu = sqrt(p/2);

%% Run models

disp("Nonnegative, NA, LOD, no missing or <LOD")
[L,S] = root_pcp_with_nan_nonnegL_LOD(X, lambda, mu, 0); 
norm(L, "Fro")
norm(S, "Fro")

disp("Nonnegative, NA, LOD, with missing")
[L,S] = root_pcp_with_nan_nonnegL_LOD(Xmissing, lambda, mu, 0); 
norm(L, "Fro")
norm(S, "Fro")

disp("Nonnegative, NA, LOD, with <LOD")
[L,S] = root_pcp_with_nan_nonnegL_LOD(Xlod, lambda, mu, 0); 
norm(L, "Fro")
norm(S, "Fro")

disp("Nonnegative, NA, LOD, with missing & <LOD")
[L,S] = root_pcp_with_nan_nonnegL_LOD(Xmissinglod, lambda, mu, 0); 
norm(L, "Fro")
norm(S, "Fro")

% disp("Nonconvex, nonnegative, NA")
% [L,S] = root_pcp_rank_r_nonnegL_with_missing(Xmissing, lambda, 1, 5); 
% norm(L, "Fro")
% norm(S, "Fro")

% disp("Nonconvex, NA")
% [Lroot,Sroot] = root_pcp_rank_r_with_missing(Xmissing, lambda, 1, 5); 
% norm(Lroot, "Fro")
% norm(Sroot, "Fro")

% disp("Nonconvex, nonnegative")
% [L1,S1] = root_pcp_rank_r_nonnegL(X, lambda, mu, 5); 
% norm(L1, "Fro")
% norm(S1, "Fro")

% disp("Nonconvex")
% [L2,S2] = root_pcp_rank_r(X, lambda, mu, 5); 
% norm(L2, "Fro")
% norm(S2, "Fro")

% disp("Nonnegative, NA")
% [Lnan_nn,Snan_nn] = root_pcp_with_nan_nonnegL(Xmissing, lambda, mu);
% norm(Lnan_nn, "Fro")
% norm(Snan_nn, "Fro")

% disp("Nonnegative, NA, LOD")
% [Lnan_nn_lod,Snan_nn_lod] = root_pcp_with_nan_nonnegL_LOD(Xmissing, lambda, mu, 0);
% norm(Lnan_nn_lod, "Fro")
% norm(Snan_nn_lod, "Fro")

% disp("Nonnegative")
% [L,S] = root_pcp_nonnegL(X, lambda, mu);
% norm(L, "Fro")
% norm(S, "Fro")

% disp("NA")
% [Lnan,Snan] = root_pcp_with_nan(Xmissing, lambda, mu);
% norm(Lnan, "Fro")
% norm(Snan, "Fro")

% disp("rootPCP")
% [L1,S1] = root_pcp(X, lambda, mu);
% norm(L1, "Fro")
% norm(S1, "Fro")


