%% Import data from text file

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 26);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21", "V22", "V23", "V24", "V25", "V26"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
lawcormat = table2array(readtable("/Users/lizzy/Principal.Component.Pursuit/root_pcp/law_cor_mat.csv", opts));
mask = table2array(readtable("/Users/lizzy/Principal.Component.Pursuit/root_pcp/law_mask.csv", opts));
mat = table2array(readtable("/Users/lizzy/Principal.Component.Pursuit/root_pcp/law_mat.csv", opts));

%% Clear temporary variables
clear opts

[n, p] = size(lawcormat);
summary = zeros(4, 8);

%%%%%%%% 1. NONCVX LOD - GOOD RUN %%%%%%%% 
lambda = 0.122;
mu = 17.6;
r = 6;

[goodL, goodS] = root_pcp_ncvx_nan_nonnegL_LOD(lawcormat, lambda, mu, r, 0);

score = norm((mat - goodL - goodS) .* mask, "fro") / norm(mat .* mask, "fro");
Lrank = rank(goodL, 1e-04);
Ssparsity = sum(sum(goodS >= -1e-04 & goodS <= 1e-04))/(n*p);
Lnorm = norm(goodL, 'fro');
Snorm = norm(goodS, 'fro');

summary(1,:) = [lambda, mu, r, score, Lrank, Ssparsity, Lnorm, Snorm];

%%%%%%% 2. NONCVX LOD - BAD RUN %%%%%%%% 
lambda = 0.272;
mu = 1.61;
r = 9;

[badL, badS] = root_pcp_ncvx_nan_nonnegL_LOD(lawcormat, lambda, mu, r, 0);

score = norm((mat - badL - badS) .* mask, "fro") / norm(mat .* mask, "fro");
Lrank = rank(badL, 1e-04);
Ssparsity = sum(sum(badS >= -1e-04 & badS <= 1e-04))/(n*p);
Lnorm = norm(badL, 'fro');
Snorm = norm(badS, 'fro');

summary(2,:) = [lambda, mu, r, score, Lrank, Ssparsity, Lnorm, Snorm];

%%%%%%%% 3. REG NONCVX - GOOD RUN %%%%%%%% 
lambda = 0.272;
mu = 13.6;
r = 4;

[regoodL, regoodS] = root_pcp_rank_r_nonnegL_with_missing(lawcormat, lambda, mu, r);

score = norm((mat - regoodL - regoodS) .* mask, "fro") / norm(mat .* mask, "fro");
Lrank = rank(regoodL, 1e-04);
Ssparsity = sum(sum(regoodS >= -1e-04 & regoodS <= 1e-04))/(n*p);
Lnorm = norm(regoodL, 'fro');
Snorm = norm(regoodS, 'fro');

summary(3,:) = [lambda, mu, r, score, Lrank, Ssparsity, Lnorm, Snorm];

%%%%%%%% 4. REG NONCVX - BAD RUN %%%%%%%% 
lambda = 0.182;
mu = 9.61;
r = 6;

[regbadL, regbadS] = root_pcp_rank_r_nonnegL_with_missing(lawcormat, lambda, mu, r);

score = norm((mat - regbadL - regbadS) .* mask, "fro") / norm(mat .* mask, "fro");
Lrank = rank(regbadL, 1e-04);
Ssparsity = sum(sum(regbadS >= -1e-04 & regbadS <= 1e-04))/(n*p);
Lnorm = norm(regbadL, 'fro');
Snorm = norm(regbadS, 'fro');

summary(4,:) = [lambda, mu, r, score, Lrank, Ssparsity, Lnorm, Snorm];

summary
