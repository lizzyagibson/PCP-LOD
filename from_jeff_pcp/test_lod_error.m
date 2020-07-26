%% Parameter definition

n = 100;
p = 20;
r = 3;      % Rank of the matrix
sigma = 0.3;  % Gaussian noise std
Delta = 0.5;  % LOD threshold, shared for all columns of the data

lambda = 1/sqrt(n);
% mu = sqrt(p/(2*log(n*p)));
mu = 1/sigma/sqrt(2*n*log(n*p));

N_try = 100;        % Number of repeats

%% Auxilliary part for setting LOD threshold Delta
% This part plots the cumulative distribution of the data entries (in D)
% generated in this way. The LOD threshold Delta can be set by referring to
% this CDF.

N = 1e6;
x_tmp = sum(rand(r,N) .* rand(r,N));
[counts,centers] = hist(x_tmp,1000);
cdf = cumsum(counts) / N;

figure
plot(centers, cdf)
xlabel('LOD threshold value')
ylabel('Proportion of data below LOD')

%% Running the algorithms

rel_err = zeros(N_try, 3);  % Matrix for storing relative error values
rel_err_belowlod = zeros(N_try, 3);  % ... for entries below LOD

for i_try = 1:N_try
    if mod(i_try,5)==0
        disp(i_try)
    end
    
    U = rand(n,r);        
    V = rand(r,p);
    L = U*V;  % Low-rank part

%     Z = randn(n,p) * sigma;  % Gaussian noise part
%     S = -(L+Z).*(L+Z<0) + (rand(n,p)<0.03) .* rand(n,p)*1;  % Sparse noise part
%     S = zeros(n,p);
%     D = L + S + Z; % Noise added
      D = L;  % No noise added
    
    D_minus1 = (D>=Delta).*D + (D<Delta)*(-1);
    D_sqrt2  = (D>=Delta).*D + (D<Delta)*(Delta/sqrt(2));
    norm_L_belowlod = norm(L.*(D<Delta), 'fro');
    norm_L = norm(L, 'fro');
    
    % Alg 1: PCA w/ LOD/sqrt(2)
    [UU,SS,VV] = svd(D_sqrt2,'econ');
    LL = UU(:,1:r) * SS(1:r,1:r) * VV(1:r,:);
    rel_err(i_try, 1) = norm(L-LL, 'fro') / norm_L;
    rel_err_belowlod(i_try, 1) = norm((L-LL).*(D<Delta), 'fro') / norm_L_belowlod;

    % Alg 2: PCP w/ LOD/sqrt(2)
    [LL,~] = pcp_lod(D_sqrt2, lambda, mu, 0);
    rel_err(i_try, 2) = norm(L-LL, 'fro') / norm_L;
    rel_err_belowlod(i_try, 2) = norm((L-LL).*(D<Delta), 'fro') / norm_L_belowlod;
    
    % Alg 3: PCP-LOD
    [LL,~] = pcp_lod(D_minus1, lambda, mu, Delta);
    rel_err(i_try, 3) = norm(L-LL, 'fro') / norm(L, 'fro');
    rel_err_belowlod(i_try, 3) = norm((L-LL).*(D<Delta), 'fro') / norm_L_belowlod;
    
end

%% Plot 1: Boxplot

figure
boxplot(rel_err, 'Notch','on', 'Labels',{'PCA w/ LOD/sqrt(2)','PCP w/ LOD/sqrt(2)','PCP-LOD'})
xlabel('Algorithms')
ylabel('Relative error in L')
title(['Delta = ',num2str(Delta)])

%% Plot 2: X-Y Scatterplot

x_scatter = rel_err(:,3)./rel_err(:,2) - 1;
y_scatter = rel_err_belowlod(:,3)./rel_err_belowlod(:,2) - 1;

figure
scatter(x_scatter, y_scatter)
line([0 0], ylim, 'Color','k');
line(xlim, [0 0], 'Color','k');
title(['Delta = ',num2str(Delta)])
xlabel('Relative error overall ratio (PCP-LOD/PCP-sqrt2 -1)')
ylabel('Relative error <LOD ratio (PCP-LOD/PCP-sqrt2 -1)')

%% Save data (if necessary)
%save(['data_',num2str(Delta),'.mat'], 'rel_err','rel_err_belowlod')
