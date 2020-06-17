%% Import data from text file
% 16-Jun-2020 12:05:32

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 10);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["PCB138", "PCB153", "PCB180", "D03", "D05", "D07", "F04", "F08", "PCBpcb", "PCBhxc"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data from wherever you have it
nhanes10 = readtable("/Users/lizzy/Principle.Component.Pursuit/Data/nhanes_10.csv", opts);
% Repeat steps with nhanes_18.csv

%% Convert to output type
nhanes10 = table2array(nhanes10);

%% Clear temporary variables
clear opts

%% Create mean vector
means10 = mean(log(nhanes10));

%% Create covariance matrix
cov10 = cov(log(nhanes10));

%% Simulate with multivariate normal
rng(1988)
mvn_sim10 = mvnrnd(means10,cov10,100);

%% exp to get multivariate log normal
mvlogn_sim10 = exp(mvn_sim10);

%% Scale multivariate log normal columns
sim10 = mvlogn_sim10 ./ std(mvlogn_sim10);
% Sim10 is the baseline matrix with no values <LOD

%% Push values <LOD to -1 if < column quantile
q50 = quantile(sim10, 0.50);
% LOD is a vector -- different LOD for each column
q50_tf = sim10 > q50;

% Create empty matrix
mix_data_lod_50 = zeros(size(sim10, 1), size(sim10, 2));

% fill matrix with sim value if >LOD and -1 if <LOD
for c = 1:size(sim10, 2)
    for r = 1:size(sim10,1)
        if q50_tf(r,c) == 1
            mix_data_lod_50(r,c) = sim10(r,c);
        else
            mix_data_lod_50(r,c) = -1;
        end
    end
end

%% mix_data_lod_50 is the simulated data w/ 50% of observations <LOD
%% Repeat with quantiles 0.10 - 0.40 to create all matrices
