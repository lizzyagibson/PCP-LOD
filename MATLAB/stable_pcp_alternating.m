function [L,S,lambda,mu] = stable_pcp_alternating(D,lambda,mu)
% D is X in demo_mixture script, lambda is 4/sqrt(m), mu is something like step size at each iteration
% lambda is how much you subtract from the singular values in the thresholding
% also 

[m,n] = size(D);
%number of rows and columns of D
S = zeros(m,n);
L = zeros(m,n);
%empty matrices of same size as D

iter = 0;
MAX_ITER = 20; 
done = false; 

while ~done,
    
    iter = iter + 1;
    
    [L,v] = singular_value_threshold( D - S, 1 / mu, [] );
	%[] parameter, 1/mu is parameter lambda
    S     = soft_thresholding( D - L, lambda / mu );
    
    obj = v + lambda * sum(sum(abs(S))) + (mu/2) * norm( D - L - S, 'fro' )^2; 
    % calculate 'fro' norm, frobenious
	% abs absolute value
	%sum(sum( gives each column sum then sum all columns
    disp( [ num2str(iter) '  Obj: ' num2str(obj) ] );
    %num2str print to output, concatenate
    if iter >= MAX_ITER,
        done = true; 
    end
end