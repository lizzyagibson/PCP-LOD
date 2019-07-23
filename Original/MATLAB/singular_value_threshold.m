function [N,v,sv] = singular_value_threshold( M, lambda, B, sv )

% B specified as [] in example
% singular value thresholding
%
% solves two problems
% 
% If nargin < 3, computes the prox operator for the nuclear norm
%   arg min_X lambda ||X||_* + .5 * ||X - M||_F^2
%
% If nargin = 3, computes a constrained prox operator
%   arg min_X lambda ||X||_* + .5 * ||X - M||_F^2   s.t.   ||X||_* <= B
%

% [m,n] = size(M); 

% if nargin < 4 || isempty(sv),
%   sv = m;
% end

% if isempty(B),
    
    %disp('SLOW');
  %  tic; 
    [U,S,V] = svd(M);
% SVD function in R
% D in R is S in MATLAB
% U is U and V is V
	
 %   toc; 
    
%     tic
%     [U,S,V] = lansvd(M,sv,'L');
%     toc
%     

   % svp = length(find(diag(S) > lambda));
   %  if svp < sv
   %     sv = min(svp + 1, m);
   % else
   %     sv = min(svp + round(0.05*m), m);
   % end    
    
    N = U * soft_thresholding(S,lambda) * V'; 
	%prime is transpose
	%S is singular values
%else
%    [U,S,V] = svd(M);
%    s = diag(S); 
%    s1 = soft_thresholding(s,lambda);
%    s2 = ProjectOntoSimplex(s,B); 
    
%    if sum(s1) < sum(s2),
%        s = s1;
%    else
%        s = s2;
%    end
    
%    N = U * diag(s) * V';
%	
end
    
v = sum(diag(soft_thresholding(S,lambda))); 
% diag converts diagonal matrix into straight vector
% takes singular value diagonal matrix and makes vector of singular values