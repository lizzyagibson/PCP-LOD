%% Parameter definition

jobid = str2num(getenv('SGE_TASK_ID'))

n = 1000;
p = 20;

sigma_vec = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
delta_vec = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
rank_vec = [1, 2, 3, 4, 5];

grid = combvec(sigma_vec, delta_vec, rank_vec)';

sigma_prop = grid(jobid, 1)
delta_prop = grid(jobid, 2)
r = grid(jobid, 3)
    
    lambda = 1/sqrt(n);
    mu = sqrt(p/(2*log(n*p)));

    N_try = 100;  % Number of repeats

    %% Running the algorithms

        rel_err = zeros(N_try, 3);           % Matrix for storing relative error values
        rel_err_belowlod = zeros(N_try, 3);  % for entries below LOD
    
    for i_try = 1:N_try
        if mod(i_try,5)==0
            disp(i_try)
        end

        U = rand(n,r);        
        V = rand(r,p);
        L = U*V;                    % Low-rank part
        
        sigma = sigma_prop * std(L(:));
        
        Z = normrnd(0,sigma,n,p);   % Gaussian noise part
        X = L + Z;                  % Noise added
        D = max(X, 0);              % Non-negative
    
        Delta = quantile(D(:), delta_prop);
        
        D_minus1 = (D>=Delta).*D + (D<Delta)*(-1);
        D_sqrt2  = (D>=Delta).*D + (D<Delta)*(Delta/sqrt(2));
        norm_L_belowlod = norm(L.*(D<Delta), 'fro');
        norm_L = norm(L, 'fro');

        % Alg 1: PCA w/ LOD/sqrt(2)
        [UU,SS,VV] = svd(D_sqrt2,'econ');
        LL = UU(:,1:r) * SS(1:r,1:r) * VV(1:r,:);
        rel_err(i_try, 1) = norm(L-LL, 'fro') / norm_L;                                % predictive error
        rel_err_belowlod(i_try, 1) = norm((L-LL).*(D<Delta), 'fro') / norm_L_belowlod; % error below lod

        % Alg 2: PCP w/ LOD/sqrt(2)
        [LL,~] = pcp_lod(D_sqrt2, lambda, mu, 0);
        rel_err(i_try, 2) = norm(L-LL, 'fro') / norm_L;
        rel_err_belowlod(i_try, 2) = norm((L-LL).*(D<Delta), 'fro') / norm_L_belowlod;

        % Alg 3: PCP-LOD
        [LL,~] = pcp_lod(D_minus1, lambda, mu, Delta);
        rel_err(i_try, 3) = norm(L-LL, 'fro') / norm(L, 'fro');
        rel_err_belowlod(i_try, 3) = norm((L-LL).*(D<Delta), 'fro') / norm_L_belowlod;
      
    end

   %% Save data 
   save(['pcp_all_lod_',num2str(delta_prop),'_rank_', num2str(r), '_sigma_', num2str(sigma_prop), '.mat'], 'rel_err','rel_err_belowlod')
