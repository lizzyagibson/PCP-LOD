function X =  proj_rank_r(Y,r)

[U,S,V] = svd(Y,'econ');

s          = diag(S);
s(r+1:end) = 0;
S_new      = diag(s);

X = U * S_new * V';
