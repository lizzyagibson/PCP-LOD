%% Data generation
n = 1000;
p = 20;
r = 3;

DeltaVec  = 0:.05:1;
numTrials = 50;
numPCPLodBetter = zeros(size(DeltaVec));
numTrialsVec = zeros(size(DeltaVec));


for i = 1:numTrials
    for j = 1:length(DeltaVec)
        
        Delta = DeltaVec(j);
        

        U = rand(n,r);        
        V = rand(r,p);
        L = U*V;  % Low-rank part

        sigma = 0;
        Z = randn(n,p) * sigma;  % Gaussian noise part

        S = zeros(n,p);
        %S = -(L+Z).*(L+Z<0) + (rand(n,p)<0.05).*rand(n,p)*0.5;  % Sparse noise part
        % S = zeros(n,p);

        D = L + S + Z;

        %% Apply algorithms

        % pcp-lod
        Din1 = (D>=Delta).*D + (D<Delta)*(-1);

        numAvail = sum( Din1 > -1, 2 );
        goodRows = numAvail > 3; 
        G = repmat(double(goodRows),1,p);

        lambda = 1/sqrt(n);
        %mu = 1/sigma/sqrt(2*n*log(n*p));
        mu = 38;
        [L1, S1] = pcp_lod(Din1, lambda, mu, Delta);
        err1 = norm(L-L1,'fro') / norm(L,'fro');
        % errBelowLOD1 = norm((D<Delta) .* (L-L1),'fro') / norm((D<Delta).* L, 'fro'); 
        errGoodRows1 = norm( G .* (L-L1),'fro') / norm( G .* L, 'fro' );
        % errGoodRows1 = norm( G .* (D<Delta) .* (L-L1),'fro') / norm( G .* (D<Delta) .* L, 'fro' ); 

        % pcp with lod/sqrt(2)
        Din2 = (D>=Delta).*D + (D<Delta)*(Delta/sqrt(2));
        [L2, S2] = pcp_lod(Din2, lambda, mu, 0);
        err2 = norm(L-L2,'fro') / norm(L,'fro');
        % errBelowLOD2 = norm((D<Delta) .* (L-L2),'fro') / norm((D<Delta).* L, 'fro'); 
        errGoodRows2 = norm( G .* (L-L2),'fro') / norm( G .* L, 'fro' ); 
        % errGoodRows2 = norm( G .* (D<Delta) .* (L-L2),'fro') / norm( G .* (D<Delta) .* L, 'fro' ); 

        % disp(' ');
        % disp(['Overall Fraction Below LOD: ' num2str(sum(sum(D < Delta)) / (n*p))]);
        % disp(' ');
        % disp(['   Relative Error in Imputed Entries, LOD-PCP  : ' num2str(errBelowLOD1)]);
        % disp(['   Relative Error in Imputed Entries, PCP root2: ' num2str(errBelowLOD2)]);
        % disp(' ');
        % disp(['Number of Missing Entries in Bad Rows : ' num2str(sum(sum( (1-G) .* double(D < Delta))))]);
        % disp(['Number of Missing Entries in Good Rows: ' num2str(sum(sum( G .* double(D < Delta))))]);
        % disp(' ');
        % disp(['   Relative Error in Good Rows, LOD-PCP   : ' num2str(errGoodRows1) ]);
        % disp(['   Relative Error in Good Rows, PCP root2 : ' num2str(errGoodRows2) ]);
        
        numTrialsVec(j) = numTrialsVec(j) + 1;
        
        if errGoodRows1 < errGoodRows2, 
            numPCPLodBetter(j) = numPCPLodBetter(j) + 1;
        end
        
        % if err1 < err2, 
        %     numPCPLodBetter(j) = numPCPLodBetter(j) + 1;
        % end
        
        figure(1);
        clf;
        plot(DeltaVec,numPCPLodBetter ./ numTrialsVec,'r-','LineWidth',3);
        pause(.1);
        
        
    end
end

