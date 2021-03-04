function [L, S] = root_pcp_rank_r_nonnegL_with_missing (D, lambda, mu, r)
% [L, S] = root_pcp_rank_r_nonnegL_with_missing (D, lambda, mu, r)
%
% Solve the following problem:
% min_{L,S}
%         ||L||_* + lambda * ||S||_1 + mu * ||P_(obs)[L+S-D]||_F
% s.t. L >= 0.
%
% This is first transformed to the problem
% min_{L1,L2,L3,S1,S2,Z}
%      ||L1||_* + lambda * ||S1||_1 + mu * ||Z||_F + I(L3>=0)
% s.t. L1 = L2
%      S1 = S2
%      Z = P_obs[ D - L2 - S2]
%      L1 = L3
% The algorithm conducts ADMM splitting as (L1,S1,Z,L3),(L2,S2).
% use nan for missing entries in the observation

[n,p] = size(D);
rho = 0.1; % Augmented Lagrangian parameter

[L1,L2,L3,S1,S2,Z,Y1,Y2,Y3,Y4] = deal(zeros(n,p));

% % % % % mask: support of observation of D
mask = ~isnan(D);
D(~mask) = 0;


MAX_ITER = 10000;
EPS_ABS = 1e-6;
EPS_REL = 1e-6;

flag_converge = 0;
% loss = zeros(MAX_ITER, 1);

% ADMM-splitting iterations
for i = 1:MAX_ITER

    % Store previous values of L2,S2
    L2_old = L2;
    L3_old = L3;
    S2_old = S2;

    % Update 1st primal variable (L1,S1,Z)
    L1 = proj_rank_r( (L2+L3-Y1/rho-Y4/rho)/2,r );
    S1 = prox_l1( S2-Y2/rho, lambda/rho );
% % % % %     Z = prox_fro( D-L2-S2-Y3/rho, mu/rho );
    Z = prox_fro( mask.*(D-L2-S2)-Y3/rho, mu/rho );

    L3 = max(L1+Y4/rho,0);
% % % % %     % Update 2nd primal variable (L2,S2)
    L2_obs = mask.*(1/3*( D-Z+2*L1-S1 + (2*Y1-Y2-Y3)/rho ));
    L2_unobs = (1-mask).*(L1+Y1/rho);
    L2 = L2_obs+L2_unobs;
    
    S2_obs = mask.*(1/3*( D-Z+2*S1-L1 + (2*Y2-Y1-Y3)/rho ));
    S2_unobs = (1-mask).*(S1+Y2/rho);
    S2 = S2_obs+S2_unobs;
    
    % Update dual variable (Y1,Y2,Y3)
    Y1 = Y1 + rho*(L1-L2);
    Y2 = Y2 + rho*(S1-S2);
% % % % %     
    Y3 = Y3 + rho*(Z-mask.*(D-L2-S2));
    Y4 = Y4 + rho*(L1-L3);
    
    %  Calculate primal & dual residuals; Update rho
% % % % %     
    res_primal = sqrt( norm(L1-L2,'fro')^2 + norm(S1-S2,'fro')^2 + ...
                       norm(Z-mask.*(D-L2-S2),'fro')^2 + norm(L1-L3,'fro')^2);
    res_dual = rho * sqrt( norm(L2+L3-L2_old-L3_old,'fro')^2 + norm(S2-S2_old,'fro')^2 + ...
                           norm(mask.*(L2-L2_old+S2-S2_old),'fro')^2 );
    if res_primal > 10 * res_dual
        rho = rho * 2;
    elseif res_dual > 10 * res_primal
        rho = rho / 2;
    end
%     figure(1);
%     heatmap(L2(1:10,1:10));
%     title(i);

    
%     % Calculate loss
%     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu*norm(L2+S2-D,'fro') ...
%         + sum(sum(Z1.*(L1-L2))) + sum(sum(Z2.*(S1-S2))) ...
%         + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((S1-S2).^2)) );
    
% % % % %     % Check stopping criteria
    thresh_primal = EPS_ABS * sqrt(4*n*p) + EPS_REL * ...
                    max([sqrt( 2*norm(L1,'fro')^2 + norm(S1,'fro')^2 + norm(Z,'fro')^2 ), ...
                         sqrt( norm(L2,'fro')^2 + norm(S2,'fro')^2 + norm(mask.*(L2+S2),'fro')^2 +norm(L3,'fro')^2), ...
                         norm(D,'fro')
                        ]);
    thresh_dual = EPS_ABS * sqrt(3*n*p) + EPS_REL * ...
                    sqrt( norm(Y1+Y4,'fro')^2 + norm(Y2,'fro')^2 + norm(Y3,'fro')^2 );
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
