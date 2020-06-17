% demo_mixtures
clear;

load('../../../Data/mixtures_data.mat');

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
[L,S, loss] = pcp_lod(X,4/sqrt(m),10, 0); 
% X is out dataset, set lambda and mu

save('../LOD_demo_output/lowrank_lod0.mat', 'L')
save('../LOD_demo_output/sparse_lod0.mat', 'S')

figure(1); imagesc(S((2417-3):(2417+5),:));

dateLabels = {'5/28/10','5/29/10','5/30/10','5/31/10','6/1/10','6/2/10','6/3/10','6/4/10','6/5/10'};
speciesLabels = {'Al','As', 'Ba', 'Bc', 'Br', 'Ca', 'Cl', 'Cr', 'Cu', 'Fe', 'K',  'Mn',  'Ni',  'Pb',  'S',  'Se',  'Si', 'Ti',  'V', 'Zn'};
saveas(gcf,'./IMAGES/fig1_lod_0.png')

set(gca,'XTick',1:size(L,2));
set(gca,'XTickLabels',speciesLabels);
set(gca,'YTickLabels',dateLabels);

figure(2);
[U,Sigma,V] = svd(L);
subplot(2,3,1);
stem(-V(:,1));
set(gca,'XTick',1:size(L,2));
set(gca,'XTickLabels',speciesLabels);
subplot(2,3,2);
stem(V(:,2));
set(gca,'XTick',1:size(L,2));
set(gca,'XTickLabels',speciesLabels);
subplot(2,3,3);
stem(V(:,3));
set(gca,'XTick',1:size(L,2));
set(gca,'XTickLabels',speciesLabels);
subplot(2,3,4);
stem(V(:,4));
set(gca,'XTick',1:size(L,2));
set(gca,'XTickLabels',speciesLabels);
subplot(2,3,5);
stem(V(:,5));
set(gca,'XTick',1:size(L,2));
set(gca,'XTickLabels',speciesLabels);
saveas(gcf,'./IMAGES/fig2_lod_0.png')