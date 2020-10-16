function [X, nuclearX] =  prox_nuclear(Y,c)

[U,S,V] = svd(Y,'econ');

S_new = sign(S) .* max(abs(S)-c, 0);
X = U * S_new * V';
nuclearX = sum(sum(abs(S_new)));
