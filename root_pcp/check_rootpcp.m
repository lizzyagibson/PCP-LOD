% demo_mixtures
clear;

addpath('/Users/lizzy/OneDrive - cumc.columbia.edu/Principal.Component.Pursuit');
load('./Data/mixtures_data.mat');

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

lambda = 1/sqrt(m); 
%weight parameter for pcp, 1 / sqrt(number of rows)

[L,S] = root_pcp_nonnegL(X, 1/sqrt(m), 10); 
[L2,S2] = root_pcp_lod(X, 1/sqrt(m), 10, 0); 
