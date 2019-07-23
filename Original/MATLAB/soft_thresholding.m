function w = soft_thresholding( v, lambda )

w = sign(v) .* max( abs(v) - lambda, 0 ); 

% .* is element wise multiplication
% sign of v, preserve the positive sign of the singular values
% shrink the value in v toward zero by lambda
% take the higher of zero of value minus lambda
