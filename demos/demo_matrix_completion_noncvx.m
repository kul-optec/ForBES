% solve a rank-constrained matrix completion problem using ForBES

close all;
clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 200; % number of rows
n = 200; % number of column of the original matrix M
d = 15000/(m*n); % density of coefficients sampled from M
r = 10; % rank of M

U = randn(m, r);
V = randn(n, r);
M = U*V';

P = sprand(m, n, d) ~= 0; % sampling pattern
B = full(M.*P);

lam = 1e0;

f = quadLoss(P(:), B(:));
g = indRankBall(m, n, 10);
x0 = zeros(m*n, 1);
opt.maxit = 1000;
opt.tol = 1e-8;
opt.Lf = 1;
opt.display = 1;

fprintf('\nFBS\n');
opt_fbs = opt;
opt_fbs.method = 'fbs';
opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, x0, [], [], opt_fbs);
fprintf('\n');
fprintf('iterations : %d\n', out_fbs.iterations);
fprintf('SVDs       : %d\n', out_fbs.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs-fpr';
out_lbfgs = forbes(f, g, x0, [], [], opt_lbfgs);
fprintf('\n');
fprintf('iterations : %d\n', out_lbfgs.iterations);
fprintf('SVDs       : %d\n', out_lbfgs.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.residual(end));
