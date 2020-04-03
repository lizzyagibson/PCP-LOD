function l = loss_lod_nn(X, D, Delta)

X_lod = (X-D) .* (D>=0) + (X-Delta) .* (D<Delta & X>Delta) + X .* (D<Delta & X<0);
l = sum(sum(X_lod.^2)) / 2;
