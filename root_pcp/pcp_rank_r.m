function [L,S,iter] = pcp_rank_r(D,gamma,r)

% solve the problem
%
%   min I_{rank(L) <=r} + gamma ||S||_1 + .5*||PO(Y-L-S)||_F^2
%
% using a proximal gradient method

[n,p] = size(D);
Om    = ~isnan(D);

D(~Om) = 0;

[L,S] = deal(zeros(n,p));
obj = inf;

t = .01;
iter = 0;

done = false; 

max_iter = 4000;

allObj = []; 

while ~done
    
    R = Om .* ( D - L - S );
    
    S_new = prox_l1( S + t * R, t*gamma );
    L_new = proj_rank_r( L + t * R, r );
    
    delta = sqrt(norm(S-S_new,'fro')^2 + norm(L-L_new,'fro')^2);
        
    obj_new = gamma * norm(S(:),1) + .5 * norm( Om .* ( D - L - S ), 'fro' )^2;
    
    if obj_new > obj, 
        t = t * .95;
    else
        t = t * 1.01;
        S = S_new;
        L = L_new;
        delta_obj = obj - obj_new;    
        obj = obj_new;   
    end
        
    iter = iter + 1;
    
    if mod(iter,100) == 0,
        disp(['Iter ' num2str(iter) '  Obj ' num2str(obj) '  Step ' num2str((delta))]);
    end
    
    if iter >= max_iter || t < 1e-10, 
        done = true;
    end
    
    allObj = [allObj; obj];
end

%figure; plot(allObj); 