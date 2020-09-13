function [L, S, flag_converge] = pcp_lod (D, lambda, mu, Delta)
%
% If the LOD threshold Delta = 0 (nothing below LOD),
% solve the following ADMM splitting problem:
% min_{L1,L2,L3,S1,S2}
%      ||L1||_* + lambda * ||S1||_1 + mu/2 * ||L2+S2-D||_F^2 + I_{L3>=0}
% s.t. L1 = L2
%      L1 = L3
%      S1 = S2.
%
% If Delta > 0, replace ||L2+S2-D||_F^2 with LOD penalty.
%
% Below-LOD data input in D should be denoted as negative values, e.g. -1.

[n,p] = size(D);
rho = 0.1; % Augmented Lagrangian parameter, init

if ~(any(size(Delta,1)==[1,n]) && any(size(Delta,2)==[1,p]))
    error('Incorrect size of Delta.')
end

[L1,L2,L3,S1,S2,Z1,Z2,Z3] = deal(zeros(n,p));

MAX_ITER = 1000;
% LOSS_THRESH = 1e-5;
% SAME_THRESH = 1e-4;

%%%%% BEGIN NEW %%%%%
EPS_ABS = 1e-4;
EPS_REL = 1e-3;

flag_converge = 0;
%%%%% END NEW %%%%%
loss = zeros(MAX_ITER, 1);

% ADMM-splitting iterations
for i = 1:MAX_ITER
    
    % Update 1st primal variable (L1,S1)
    [L1, nuclearL1] = prox_nuclear( (L2+L3-(Z1+Z2)/rho)/2, 1/2/rho );
    S1 = prox_l1( S2-Z3/rho, lambda/rho );
    
    % Update 2nd primal variable (L2,L3,S2)
    %%%%% BEGIN NEW %%%%%
    L2_old = L2;
    L3_old = L3;
    S2_old = S2;
    %%%%% END NEW %%%%%
    
    L2_opt1 = (mu*rho*D + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho+rho^2);
    L2_opt2 = L1 + Z1/rho;
    L2_opt3 = (mu*rho*Delta + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho+rho^2);
    L2_opt4 = ((mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho+rho^2);
    L2_new = L2_opt1 .* (D>=0) + L2_opt2 .* (D<0 & L2+S2>=0 & L2+S2<=Delta) + ...
        L2_opt3 .* (D<0 & L2+S2>Delta) + L2_opt4 .* (D<0 & L2+S2<0);
    S2_opt1 = (mu*rho*D + (mu+rho)*Z3 - mu*Z1 + (mu+rho)*rho*S1 - mu*rho*L1) / (2*mu*rho+rho^2);
    S2_opt2 = S1 + Z3/rho;
    S2_opt3 = (mu*rho*Delta + (mu+rho)*Z3 - mu*Z1 + (mu+rho)*rho*S1 - mu*rho*L1) / (2*mu*rho+rho^2);
    S2_opt4 = ((mu+rho)*Z3 - mu*Z1 + (mu+rho)*rho*S1 - mu*rho*L1) / (2*mu*rho+rho^2);
    S2 = S2_opt1 .* (D>=0) + S2_opt2 .* (D<0 & L2+S2>=0 & L2+S2<=Delta) + ...
        S2_opt3 .* (D<0 & L2+S2>Delta) + S2_opt4 .* (D<0 & L2+S2<0);
    L2 = L2_new;
    L3 = max(L1+Z2/rho, 0);
    
    %%%%% BEGIN NEW %%%%%
    % Calculate primal & dual residuals; Update rho
    res_primal = sqrt( norm(L1-L2,'fro')^2 + norm(L1-L3,'fro')^2 + norm(S1-S2,'fro')^2 );
    res_dual = rho * sqrt( norm(L2+L3-L2_old-L3_old,'fro')^2 + norm(S2-S2_old,'fro')^2 );
    if res_primal > 10 * res_dual
        rho = rho * 2;
    elseif res_dual > 10 * res_primal
        rho = rho / 2;
    end 
    %%%%% END NEW %%%%%
    
    % Update dual variable (Z1,Z2,Z3)
    Z1 = Z1 + rho*(L1-L2);
    Z2 = Z2 + rho*(L1-L3);
    Z3 = Z3 + rho*(S1-S2);

    loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu*loss_lod(L2+S2,D,Delta) ...
        + sum(sum(Z1.*(L1-L2))) + sum(sum(Z2.*(L1-L3))) + sum(sum(Z3.*(S1-S2))) ...
        + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((L1-L3).^2)) + sum(sum((S1-S2).^2)) );
    
    %%%%% BEGIN NEW %%%%%
    % Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(3*n*p) + EPS_REL * ...
                    max([sqrt( norm(L1,'fro')^2 * 2 + norm(S1,'fro')^2 ), ...
                         sqrt( norm(L2,'fro')^2 + norm(L3,'fro')^2 + norm(S2,'fro')^2 )
                        ]);
    thresh_dual = EPS_ABS * sqrt(2*n*p) + EPS_REL * ...
                    sqrt( norm(Z1+Z2,'fro')^2 + norm(Z3,'fro')^2 );
    if res_primal < thresh_primal && res_dual < thresh_dual
        flag_converge = 1;
        disp(['Converged in ',num2str(i),' iterations.']);
        break
    end
    %%%%% END NEW %%%%%

end

L = L3; % (L1+L2+L3) / 3;
S = S1 % (S1+S2) / 2;

%%%%% BEGIN NEW %%%%%
if flag_converge == 0
    disp('Did not converge.');
end
%%%%% END NEW %%%%%

% figure
% plot(loss(1:i))

return

