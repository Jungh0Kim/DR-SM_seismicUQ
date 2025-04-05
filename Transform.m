% The variable x_S can be visualized and used directly, but not recommended.

damp.mu = 0.03;       damp.cov = 0.2;
E.mu = 200000;        E.cov = 0.05;
fy_b.mu = 248;        fy_b.cov = 0.1;
fy_c.mu = 345;        fy_c.cov = 0.1;
Ha_b.mu = 0.01;       Ha_b.cov = 0.2;
Ha_c.mu = 0.01;       Ha_c.cov = 0.2;

[damp.lambda, damp.zeta] = Logn_para(damp.mu, damp.mu*damp.cov);
[E.lambda, E.zeta] = Logn_para(E.mu, E.mu*E.cov);
[fy_b.lambda, fy_b.zeta] = Logn_para(fy_b.mu, fy_b.mu*fy_b.cov);
[fy_c.lambda, fy_c.zeta] = Logn_para(fy_c.mu, fy_c.mu*fy_c.cov);
[Ha_b.lambda, Ha_b.zeta] = Logn_para(Ha_b.mu, Ha_b.mu*Ha_b.cov);
[Ha_c.lambda, Ha_c.zeta] = Logn_para(Ha_c.mu, Ha_c.mu*Ha_c.cov);

x_S_train2 = zeros(n_train, size(x_S_train,2));
x_S_train2(:,1) = logninv(normcdf(x_S_train(:,1)), damp.lambda, damp.zeta);
x_S_train2(:,2) = logninv(normcdf(x_S_train(:,2)), E.lambda, E.zeta);
x_S_train2(:,3) = logninv(normcdf(x_S_train(:,3)), fy_b.lambda, fy_b.zeta);
x_S_train2(:,4) = logninv(normcdf(x_S_train(:,4)), fy_c.lambda, fy_c.zeta);
x_S_train2(:,5) = logninv(normcdf(x_S_train(:,5)), Ha_b.lambda, Ha_b.zeta);
x_S_train2(:,6) = logninv(normcdf(x_S_train(:,6)), Ha_c.lambda, Ha_c.zeta);
x_S_train = x_S_train2;

x_S_test2 = zeros(n_test, size(x_S_test,2));
x_S_test2(:,1) = logninv(normcdf(x_S_test(:,1)), damp.lambda, damp.zeta);
x_S_test2(:,2) = logninv(normcdf(x_S_test(:,2)), E.lambda, E.zeta);
x_S_test2(:,3) = logninv(normcdf(x_S_test(:,3)), fy_b.lambda, fy_b.zeta);
x_S_test2(:,4) = logninv(normcdf(x_S_test(:,4)), fy_c.lambda, fy_c.zeta);
x_S_test2(:,5) = logninv(normcdf(x_S_test(:,5)), Ha_b.lambda, Ha_b.zeta);
x_S_test2(:,6) = logninv(normcdf(x_S_test(:,6)), Ha_c.lambda, Ha_c.zeta);
x_S_test = x_S_test2;

% figure()
% histogram(x_S_train(:,2))