function [L, S] = root_pcp_nonnegL (D, lambda, mu, Delta)
% [L, S] = sqrt_pcp_LOD( D, lambda, mu, Delta)
%
% If the LOD threshold Delta = 0 (nothing below LOD),
% solve the following ADMM splitting problem:
% 
% min_{L,S}
%         ||L||_* + lambda * ||S||_1 + mu * ||L+S-D||_F
% s.t. L >= 0.
%
% This is first transformed to the problem:
%
% min_{L1,L2,L3,S1,S2,Z}
%      ||L1||_* + lambda * ||S1||_1 + mu * ||Z||_F + I(L3>=0)
%
% s.t. L1 = L2
%      S1 = S2
%      L2 + S2 + Z = D.
%      L3 = L2

% If Delta > 0, replace ||L2+S2-D||_F with LOD penalty.
%    If value < 0, replace ||L2+S2-D||_F with 0.
%    If value > LOD, replace ||L2+S2-D||_F with ||L2+S2-Delta||_F.
%
% Below-LOD data input in D should be denoted as negative values, e.g. -1.

% The algorithm conducts ADMM splitting as (L1,S1,Z,L3),(L2,S2).

[n,p] = size(D);
rho = 0.1; % Augmented Lagrangian parameter

if ~(any(size(Delta,1)==[1,n]) && any(size(Delta,2)==[1,p]))
    error('Incorrect size of Delta.')
end

[L1,L2,L3,S1,S2,Z,Y1,Y2,Y3,Y4] = deal(zeros(n,p));

MAX_ITER = 10000;
EPS_ABS = 1e-6;
EPS_REL = 1e-6;

flag_converge = 0;
% loss = zeros(MAX_ITER, 1);

% ADMM-splitting iterations
for i = 1:MAX_ITER
    
    % Store previous values of L2,S2
    L2_old = L2;
    S2_old = S2;
    
    % Update 1st primal variable (L1,S1,Z,L3)
    [L1, ~] = prox_nuclear( L2-Y1/rho, 1/rho );
    S1 = prox_l1( S2-Y2/rho, lambda/rho );
    Z = prox_fro( D-L2-S2-Y3/rho, mu/rho );
    L3 = max(L2-Y4/rho, 0);
    
    % Update 2nd primal variable (L2,S2)
    % This step depends on LOD
    % 4 different derivatives
    term1 = L1+L3+D-Z + (Y1+Y4-Y3)/rho;
    term2 = S1+D-Z + (Y2-Y3)/rho;
    
    term1_LOD = L1+L3+Delta-Z + (Y1+Y4-Y3)/rho;
    term2_LOD = S1+Delta-Z + (Y2-Y3)/rho;
    
    term1_0 = L1+L3+0-Z + (Y1+Y4-Y3)/rho;
    term2_0 = S1+0-Z + (Y2-Y3)/rho;
    
    L2_opt1 = (2*term1 - term2) / 5;
    S2_opt1 = (-term1 + 3*term2) / 5;
    
    L2_opt2 = L1 + Y1/rho; % error term drops out
    S2_opt2 = S1 + Y2/rho;
    
    L2_opt3 = (2*term1_LOD - term2_LOD) / 5;
    S2_opt3 = (-term1_LOD + 3*term2_LOD) / 5;
    
    L2_opt4 = (2*term1_0 - term2_0) / 5;
    S2_opt4 = (-term1_0 + 3*term2_0) / 5;
    
    % Keep the derivative step that agrees with the data
    L2 = L2_opt1 .* (D>=0) + L2_opt2 .* (D<0 & L2+S2>=0 & L2+S2<=Delta) + ...
        L2_opt3 .* (D<0 & L2+S2>Delta) + L2_opt4 .* (D<0 & L2+S2<0);
    S2 = S2_opt1 .* (D>=0) + S2_opt2 .* (D<0 & L2+S2>=0 & L2+S2<=Delta) + ...
        S2_opt3 .* (D<0 & L2+S2>Delta) + S2_opt4 .* (D<0 & L2+S2<0);
    
    % Update dual variable (Y1,Y2,Y3)
    Y1 = Y1 + rho*(L1-L2);
    Y2 = Y2 + rho*(S1-S2);
    Y3 = Y3 + rho*(L2+S2+Z-D);
    Y4 = Y4 + rho*(L3-L2);
    
    %  Calculate primal & dual residuals; Update rho
    res_primal = sqrt( norm(L1-L2,'fro')^2 + norm(S1-S2,'fro')^2 + ...
                       norm(Z+L2+S2-D,'fro')^2 + norm(L3-L2,'fro')^2 );
    res_dual = rho * sqrt( 2 * norm(L2-L2_old,'fro')^2 + norm(S2-S2_old,'fro')^2 + ...
                           norm(L2-L2_old+S2-S2_old,'fro')^2 );
    if res_primal > 10 * res_dual
        rho = rho * 2;
    elseif res_dual > 10 * res_primal
        rho = rho / 2;
    end 
    
%     % Calculate loss
%     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu*norm(L2+S2-D,'fro') ...
%         + sum(sum(Y1.*(L1-L2))) + sum(sum(Y2.*(S1-S2))) + sum(sum(Y3.*(L2+S2+Z-D))) + sum(sum(Y4.*(L3-L2))) ...
%         + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((S1-S2).^2)) + sum(sum(L2+S2+Z-D)) + sum(sum(L3-L2)) );
    
    % Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(4*n*p) + EPS_REL * ...
                    max([sqrt( norm(L1,'fro')^2 + norm(S1,'fro')^2 + norm(Z,'fro')^2 ) + norm(L3,'fro')^2, ...
                         sqrt( 2 * norm(L2,'fro')^2 + norm(S2,'fro')^2 + norm(L2+S2,'fro')^2 ), ...
                         norm(D,'fro')
                        ]);
    thresh_dual = EPS_ABS * sqrt(4*n*p) + EPS_REL * ...
                    sqrt( norm(Y1,'fro')^2 + norm(Y2,'fro')^2 + norm(Y3,'fro')^2 + norm(Y4,'fro')^2 );
    if res_primal < thresh_primal && res_dual < thresh_dual
        flag_converge = 1;
        disp(['Converged in ',num2str(i),' iterations.']);
        break
    end
end

L = (L1+L2+L3) / 3;
S = (S1+S2) / 2;
if flag_converge == 0
    disp('Did not converge.');
end

% figure
% plot(loss(1:i))

return
