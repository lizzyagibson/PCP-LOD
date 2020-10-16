function Y = prox_fro(X,c)

n = norm(X,'fro');

if n <= c
	Y = zeros(size(X));
else
	Y = (1 - c/n) * X;
end
