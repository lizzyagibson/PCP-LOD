function l = loss_lod(X, D, Delta)

X_lod = (X-D) .* (D>=0) + (X-Delta) .* (D<0 & X>Delta) + X .* (D<0 & X<0);
l = sum(sum(X_lod.^2)) / 2;
