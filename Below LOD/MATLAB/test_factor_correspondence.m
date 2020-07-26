% demo factor correspondence

m = 30;
n = 5;

B = randn(m,n);

% normalize the columns to have unit L2 norm
B = B ./ repmat(sqrt(sum(B.*B,1)),m,1); 

I = eye(n);
p = randperm(n);
P = I(:,p);

A = B * P;

[e,P_hat] = factor_correspondence(A,B,true);

disp(['Error: ' num2str(e)]);
figure(1);
clf;
subplot(1,2,1);
imagesc(P);
subplot(1,2,2);
imagesc(P_hat);


pause; 

% now do it again with noise
A = B * P;
A = A + .25 * randn(size(A));
A = A ./ repmat(sqrt(sum(A.*A,1)),m,1); 

[e,P_hat] = factor_correspondence(A,B,true);

disp(['Error: ' num2str(e)]);
figure(1);
clf;
subplot(1,2,1);
imagesc(P);
subplot(1,2,2);
imagesc(P_hat);