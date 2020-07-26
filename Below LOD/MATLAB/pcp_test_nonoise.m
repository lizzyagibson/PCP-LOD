%% Data generation
r = getenv('SGE_TASK_ID')

n = 1000;
p = 20;

% Loop these:
% Rank of the matrix
Delta_prop = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]; 
% LOD threshold, shared for all columns of the data

numTrials = 100;
numPCPLodBetter = zeros(size(Delta_prop));
numTrialsVec = zeros(size(Delta_prop));

rel_err = zeros(numTrials, 2);           % Matrix for storing relative error values
rel_err_belowlod = zeros(numTrials, 2);  % for entries below LOD
rel_err_safe = zeros(numTrials, 2);  

for i = 1:numTrials
    for j = 1:length(Delta_prop)
        
        prop = Delta_prop(j);
        
        U = rand(n,r);        
        V = rand(r,p);
        L = U*V;  % Low-rank part

        sigma = 0;
        Z = randn(n,p) * sigma;  % Gaussian noise part

        S = zeros(n,p);

        D = L + S + Z;

        %% Apply algorithms
        Delta = quantile(D(:), prop);

        norm_L = norm(L, 'fro');

        % pcp-lod
        Din1 = (D>=Delta).*D + (D<Delta)*(-1);

        numAvail = sum( Din1 > -1, 2 );
        goodRows = numAvail > r; 
        G = repmat(double(goodRows),1,p);

        lambda = 1/sqrt(n);
        %mu = 1/sigma/sqrt(2*n*log(n*p));
        mu = sqrt(p/(2*log(n*p)));

        [L1, S1] = pcp_lod(Din1, lambda, mu, Delta);
        err1 = norm(L-L1,'fro') / norm(L,'fro');
        % errGoodRows1 = norm( G .* (D<Delta) .* (L-L1),'fro') / norm( G .* (D<Delta) .* L, 'fro' ); 

        rel_err(i, 1) = norm(L-L1, 'fro') / norm_L;
        rel_err_belowlod(i, 1) = norm((D<Delta) .* (L-L1),'fro') / norm((D<Delta).* L, 'fro'); 
        rel_err_safe(i, 1) = norm( G .* (L-L1),'fro') / norm( G .* L, 'fro' );

        % pcp with lod/sqrt(2)
        Din2 = (D>=Delta).*D + (D<Delta)*(Delta/sqrt(2));
        [L2, S2] = pcp_lod(Din2, lambda, mu, 0);
        err2 = norm(L-L2,'fro') / norm(L,'fro');
        % errGoodRows2 = norm( G .* (D<Delta) .* (L-L2),'fro') / norm( G .* (D<Delta) .* L, 'fro' ); 
        
        rel_err(i, 2) = norm(L-L2, 'fro') / norm_L;
        rel_err_belowlod(i, 2) = norm((D<Delta) .* (L-L2),'fro') / norm((D<Delta).* L, 'fro'); 
        rel_err_safe(i, 2) = norm( G .* (L-L2),'fro') / norm( G .* L, 'fro' ); 
        
    end
    %% Save data (if necessary)
    save(['no_noise_lod_',num2str(delta_prop),'_rank_', num2str(r),  '.mat'], 'rel_err','rel_err_belowlod', 'rel_err_safe')

end

