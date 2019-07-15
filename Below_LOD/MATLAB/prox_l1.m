function X = prox_l1(Y,c)

X = sign(Y) .* max(abs(Y)-c, 0);
