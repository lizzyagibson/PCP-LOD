function [L, S] = pcp_lod (D, lambda, mu, Delta)
%
% If the LOD threshold Delta = 0, solve the following ADMM splitting problem:
% min_{L1,L2,L3,S1,S2}
%      ||L1||_* + lambda * ||S1||_1 + mu/2 * ||L2+S2-D||_F^2 + I_{L3>=0}
% s.t. L1 = L2
%      L1 = L3
%      S1 = S2.
%
% If Delta is not 0, replace ||L2+S2-D||_F^2 with LOD penalty.
%
% Below-LOD data input in D should be denoted as negative values, e.g. -1.

% % /Test/
% Delta = 0.02;
% U = rand(5,1);
% V = rand(1,5);
% D = U*V .* (rand(5)>0.05);
% lambda = 0.5;
% mu = 20;
% % /Test/

[m,n] = size(D);
rho = 1; % Augmented Lagrangian coefficient

[L1,L2,L3,S1,S2,Z1,Z2,Z3] = deal(zeros(m,n));

MAX_ITER = 100;
LOSS_THRESH = 1e-5;
SAME_THRESH = 1e-4;
loss = zeros(MAX_ITER, 1);

for i = 1:MAX_ITER
    [L1, nuclearL1] = prox_nuclear( (L2+L3-(Z1+Z2)/rho)/2, 1/2/rho );
    S1 = prox_l1( S2-Z3/rho, lambda/rho );
    
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
    % The code block above takes LOD into account.
    % The code block commented out below does not take LOD into account
%     L2 = (mu*rho*D + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho+rho^2);
%     S2 = (mu*rho*D + (mu+rho)*Z3 - mu*Z1 + (mu+rho)*rho*S1 - mu*rho*L1) / (2*mu*rho+rho^2);
    L3 = max(L1+Z2/rho, 0);
    
    Z1 = Z1 + rho*(L1-L2);
    Z2 = Z2 + rho*(L1-L3);
    Z3 = Z3 + rho*(S1-S2);

    loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu*loss_lod(L2+S2,D,Delta) ...
        + sum(sum(Z1.*(L1-L2))) + sum(sum(Z2.*(L1-L3))) + sum(sum(Z3.*(S1-S2))) ...
        + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((L1-L3).^2)) + sum(sum((S1-S2).^2)) );
    % The code block above takes LOD into account.
    % The code block commented out below does not take LOD into account
%     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu/2*sum(sum((L2+S2-D).^2)) ...
%         + sum(sum(Z1.*(L1-L2))) + sum(sum(Z2.*(L1-L3))) + sum(sum(Z3.*(S1-S2))) ...
%         + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((L1-L3).^2)) + sum(sum((S1-S2).^2)) );
    
    if i ~= 1 && abs(loss(i-1)-loss(i)) < LOSS_THRESH ...
            && is_same(SAME_THRESH, L1,L2,L3) && is_same(SAME_THRESH, S1,S2)
        break
    end
end

L = (L1+L2+L3) / 3;
S = (S1+S2) / 2;

% figure
% plot(loss(1:i))
