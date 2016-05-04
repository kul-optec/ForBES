% solve a sparse robust regression problem using ForBES

close all;
clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 5000;
n = 20000;
x_orig = sprandn(n, 1, 100/n);
A = sprandn(m, n, 300/n);
outl = rand(m,1) < 0.01;
% add small noise for everyone + large noise for outliers
b = A*x_orig + randn(m, 1)/10 + outl.*randn(m, 1)*10;
lam_max = norm(A'*b,'inf');
lam = 0.1*lam_max;
% since we know what small/large noise means, would do cross-validation otherwise I guess
del = 1;

f = huberLoss(del);
aff = {A, -b};
g = l1Norm(lam);
x0 = zeros(n, 1);
opt.maxit = 10000;
opt.tol = 1e-8;
opt.display = 1;

fprintf('\nFast FBS\n');
opt_fbs = opt;
opt_fbs.solver = 'fbs';
opt_fbs.variant = 'fast';
out = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('matvecs    : %d\n', out.operations.C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('matvecs    : %d\n', out.operations.C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));


