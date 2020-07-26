function [e,Pi] = factor_correspondence(A,B,nn)

% factor_correspondence
%
%    find a matching between two sets of factors (presented as columns of
%    matrices $A$ and $B$) that minimizes the sum of the squared errors. 
%
%    Assumes that $A$ and $B$ have columns of unit $L2$ norm. Looks for a
%    permutation matrix (or signed permutation matrix) $\Pi$ such that 
%
%       || A - B * Pi ||_F^2 
%
%    is minimized.
%
%    Uses an (exact) Lp relaxation, solved via the cvx package.
%
%    Inputs
%       A, B matrices of the same size, with unit L2 norm columns
%       nn   boolean, optional, default true. If nn, we search for a
%       *permutation* matrix Pi. If nn is false, we allow *signed
%       permutations* which can multiply the factors by -1. 
%
%   Outputs
%       e -- the sum of squared errors under the best calibration \Pi
%       Pi -- the permutation or signed permutation matrix
%
%    June 2020 John Wright, jw2966@columbia.edu
%
%%%%


if nargin < 3, 
    nn = true; 
end

if ~all(size(A) == size(B))
    error('current code only works when A and B are the same size');
end

G = B' * A; 
[n,~] = size(G);

if nn

    cvx_begin
        variable Pi(n,n)

        maximize sum(sum(Pi .* G))
        subject to 

            Pi >= 0;
            sum(Pi,1) == 1;
            sum(Pi,2) == 1;

    cvx_end
    
else
    % allow sign flips 
    
    cvx_begin 
        variable Pi(n,n) 
        
        maximize sum(sum(Pi .* G))
        subject to 
            norms(Pi,1,1) <= 1;
            norms(Pi,1,2) <= 1;
    cvx_end
end

e = norm(B,'fro')^2 + norm(A,'fro')^2 - 2 * cvx_optval;
