% demo_mixtures
clear;

% Initialize variables.
filename = '/Users/lizzy/Principle.Component.Pursuit/Below_LOD/R/BLOD_airpol_data/mix_data_lod_20.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
Al = dataArray{:, 1};
As = dataArray{:, 2};
Ba = dataArray{:, 3};
bc = dataArray{:, 4};
Br = dataArray{:, 5};
Ca = dataArray{:, 6};
Cl = dataArray{:, 7};
Cr = dataArray{:, 8};
Cu = dataArray{:, 9};
Fe = dataArray{:, 10};
K = dataArray{:, 11};
Mn = dataArray{:, 12};
Ni = dataArray{:, 13};
Pb = dataArray{:, 14};
S = dataArray{:, 15};
Se = dataArray{:, 16};
Si = dataArray{:, 17};
Ti = dataArray{:, 18};
V = dataArray{:, 19};
Zn = dataArray{:, 20};

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%% BEGIN PCP
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta = vector of 20% quantiles for 20 variables
delta20 = [0.0268, 0, 0, 0.3461307, 0, 0.0182, 0.0011, 0, 0.0011, ...
    0.0348, 0.0211, 0, 5e-04, 0.0027, 0.41254, 0, 0.03334, 0.0017, ...
    8e-04, 0.0046]

%[L,S] = constrained_completion(X,Omega,2000,1e-8,1e-8,1,lambda,lambda,0,0,0,0,0,0); 
[L,S] = pcp_lod(X, 4/sqrt(m), 10, delta20); 
% X is out dataset, set lambda and mu

save('../LOD_demo_output/lowrank_lod20.mat', 'L')
save('../LOD_demo_output/sparse_lod20.mat', 'S')

figure(1); imagesc(S((2417-3):(2417+5),:));

dateLabels = {'5/28/10','5/29/10','5/30/10','5/31/10','6/1/10','6/2/10','6/3/10','6/4/10','6/5/10'};
speciesLabels = {'Al','As', 'Ba', 'Bc', 'Br', 'Ca', 'Cl', 'Cr', 'Cu', 'Fe', 'K',  'Mn',  'Ni',  'Pb',  'S',  'Se',  'Si', 'Ti',  'V', 'Zn'};
saveas(gcf,'../../IMAGES/fig1_lod_20.png')
% gcf means current figure

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
saveas(gcf,'../../IMAGES/fig2_lod_20.png')
