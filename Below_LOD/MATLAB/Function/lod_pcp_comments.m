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
%% Below-LOD data input in D should be denoted as negative values, e.g. -1.

% % /Test/
% Delta = 0.02;
% U = rand(5,1);
% V = rand(1,5);
% D = U*V .* (rand(5)>0.05);
% lambda = 0.5;
% mu = 20;
% % /Test/

[m,n] = size(D);
rho = 1; % Augmented Lagrangian coefficient aka rate

[L1,L2,L3,S1,S2,Z1,Z2,Z3] = deal(zeros(m,n));

MAX_ITER = 100;
% Should this be rank / dimension / condition number dependent?
LOSS_THRESH = 1e-5;
SAME_THRESH = 1e-4;
% Thresholds
loss = zeros(MAX_ITER, 1);
% could we have an autodifferentiation version of this code?
for i = 1:MAX_ITER
    [L1, nuclearL1] = prox_nuclear( (L2 + L3 - (Z1+Z2)/rho )/2, 1/2/rho );
    % L, Z, S all start at zero, and change each iteration
    % Prox_nuc is singular value thresholding
    % L is low rank matrix
    S1 = prox_l1( S2-Z3/rho, lambda/rho );
    % S is sparse matrix
    L2_opt1 = (mu*rho*D     + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2);
    L2_opt2 = L1 + Z1/rho;
    L2_opt3 = (mu*rho*Delta + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2);
    L2_opt4 = (               (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho + rho^2);
    L2_new =  L2_opt1 .* (D>=0)  ...
            + L2_opt2 .* (D<0 & L2 + S2>=0 & L2 + S2<=Delta) ...
            + L2_opt3 .* (D<0 & L2 + S2>Delta)  ...
            + L2_opt4 .* (D<0 & L2 + S2<0);
    % For positive data, use opt 1
    % L2_new becomes whichever of the 4 meets the conditions
    S2_opt1 = (mu*rho*D     + (mu+rho)*Z3 - mu*Z1 + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2);
    S2_opt2 = S1 + Z3/rho;
    S2_opt3 = (mu*rho*Delta + (mu+rho)*Z3 - mu*Z1 + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2);
    S2_opt4 = (               (mu+rho)*Z3 - mu*Z1 + (mu + rho)*rho*S1 - mu*rho*L1) / (2*mu*rho + rho^2);
    S2 =      S2_opt1 .* (D>=0) ...
            + S2_opt2 .* (D<0 & L2 + S2>=0 & L2+S2<=Delta) ...
            + S2_opt3 .* (D<0 & L2 + S2>Delta) ...
            + S2_opt4 .* (D<0 & L2 + S2<0);
    % For positive data, use opt 1
    % S2 becomes whichever of the 4 meets the conditions
    L2 = L2_new;
    % The code block above takes LOD into account.
    % The code block commented out below does not take LOD into account
%     L2 = (mu*rho*D + (mu+rho)*Z1 - mu*Z3 + (mu+rho)*rho*L1 - mu*rho*S1) / (2*mu*rho+rho^2);
%     S2 = (mu*rho*D + (mu+rho)*Z3 - mu*Z1 + (mu+rho)*rho*S1 - mu*rho*L1) / (2*mu*rho+rho^2);
    L3 = max(L1+Z2/rho, 0);
    % Initially zero at first iteration
    Z1 = Z1 + rho*(L1-L2);
    Z2 = Z2 + rho*(L1-L3);
    Z3 = Z3 + rho*(S1-S2);
    % Z accumulate differnces between L and L and between S and S

    loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) 
    + mu*loss_lod(L2+S2,D,Delta) ...
        + sum(sum(Z1.*(L1-L2))) ...
        + sum(sum(Z2.*(L1-L3))) ...
        + sum(sum(Z3.*(S1-S2))) ...
        + rho/2 * ( sum(sum((L1-L2).^2))  ...
        + sum(sum((L1-L3).^2))  ...
        + sum(sum((S1-S2).^2)) );
    % The code block above takes LOD into account.
    % The code block commented out below does not take LOD into account
%     loss(i) = nuclearL1 + lambda*sum(sum(abs(S1))) + mu/2*sum(sum((L2+S2-D).^2)) ...
%         + sum(sum(Z1.*(L1-L2))) + sum(sum(Z2.*(L1-L3))) + sum(sum(Z3.*(S1-S2))) ...
%         + rho/2 * ( sum(sum((L1-L2).^2)) + sum(sum((L1-L3).^2)) + sum(sum((S1-S2).^2)) );
    
    if i ~= 1 && abs(loss(i-1)-loss(i)) < LOSS_THRESH ...
            && is_same(SAME_THRESH, L1,L2,L3) && is_same(SAME_THRESH, S1,S2)
            % Convergence criteria!
        break
    end
end

L = (L1+L2+L3) / 3;
S = (S1+S2) / 2;

% figure
% plot(loss(1:i))

%%%%%%
%%%%%%
%%%%%%

function flag = is_same(SAME_THRESH, varargin)
% Is the difference in a sequence of matrices is above some noise threshold
% Then FALSE
% Each max() shrinks the rank down by one
% Comparing the inputs
flag = true;
for i = 1 : nargin-1
% nargin returns the number of function input arguments given in the call to 
% the currently executing function.
    for j = i+1 : nargin-1
        if max(max(abs(varargin{i}-varargin{j}))) > SAME_THRESH
            flag = false;
            return
        end
    end
end

%%%%%%
%%%%%%
%%%%%%

function l = loss_lod(X, D, Delta)
% D is the original data
% X is the new thing
X_lod = (X-D) .* (D>=0) 
% Pointwise boolean operation tricks
% D>=0 will spit out 0/1 (no/yes)
% If D_ij >= 0, then it = (X - D)_ij, else zero
      + (X-Delta) .* (D<0 & X>Delta)
% If D_ij < 0 AND X_ij > Delta, then X_ij - Delta, else zero
% D CANNOT BE < 0
      + X .* (D<0 & X<0);
% If D_ij < 0 AND X_ij < 0, then X, else zero
% D CANNOT BE < 0

l = sum(sum(X_lod.^2)) / 2;
% L2 norm squared

% Any D_ij < 0 AND X_ij < Delta AND > 0 are treated as equal
% Want to minimize the loss
% Minimize discrepancy for valid data
% Want to shrink negative things

%%%%%%
%%%%%%
%%%%%%

function X = prox_l1(Y,c)
% This is only the threshold
% If it is smaller than c, you get zero
% prox_l1 is soft thresholding
X = sign(Y) .* max(abs(Y)-c, 0);

%%%%%%
%%%%%%
%%%%%%

function [X, nuclearX] =  prox_nuclear(Y,c)
% Nuclear norm is the L1 norm of the singular values
% This encourages the matrix to be low rank bc a sparse singular vector 
% is the def of low rank
% aka it is mostly zeroes
[U,S,V] = svd(Y);

S_new = sign(S) .* max(abs(S)-c, 0);
% Threshold the singular values
% If the singular value is less than c, push it to zero
X = U * S_new * V';
% X is the truncated of the original
% Multiply the thresholded SVD components back together
nuclearX = sum(sum(abs(S_new)));
% This is the L1 of the truncated singular values

%%%%%%
%%%%%%
%%%%%%
