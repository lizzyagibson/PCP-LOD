% demo_mixtures
clear;

load('mixtures_data.mat');

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
%X = normalize_columns(X); 
% subtract mean and divide by std dev

%Xcentered = X; 
%colMeans  = sum(X,1) / m; 
%Xcentered = X - repmat(colMeans,m,1);
%colStd = sqrt( sum( Xcentered.^2, 1 ) / m );
%X = Xcentered ./ repmat(colStd,m,1); 

Omega = ~isnan(X); 
%not used right now

lambda = 1/sqrt(m); 
%weight parameter for pcp, 1 / sqrt(number of rows)

%[L,S] = constrained_completion(X,Omega,2000,1e-8,1e-8,1,lambda,lambda,0,0,0,0,0,0); 
[L,S, loss] = pcp_lod(X, 4/sqrt(m), 10, 0); 
% X is out dataset, set lambda and mu

save('lowrank_test.mat', 'L')
save('sparse_test.mat', 'S')
