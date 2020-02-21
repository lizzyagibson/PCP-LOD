%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: /Users/lizzy/Principle.Component.Pursuit/Below_LOD/R/BLOD_airpol_data/mix_data_lod_10_matrixlod.csv
%
% Auto-generated by MATLAB on 21-Feb-2020 13:31:39

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 20);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Al", "As", "Ba", "bc", "Br", "Ca", "Cl", "Cr", "Cu", "Fe", "K", "Mn", "Ni", "Pb", "S", "Se", "Si", "Ti", "V", "Zn"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
mixdatalod10matrixlod = readtable("/Users/lizzy/Principle.Component.Pursuit/Below_LOD/R/BLOD_airpol_data/mix_data_lod_10_matrixlod.csv", opts);


%% Clear temporary variables
clear opts