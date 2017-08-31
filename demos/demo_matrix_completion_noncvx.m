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

f = forbes.functions.SqrNormL2(P);
aff = {1, -B};
g = forbes.functions.IndRankBall(5);
x0 = zeros(m, n);
opt.maxit = 1000;
opt.tol = 1e-6;
opt.Lf = 1;
opt.display = 1;

fprintf('\nFBS\n');
opt_fbs = opt;
opt_fbs.solver = 'fbs';
opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('\n');
fprintf('iterations : %d\n', out_fbs.solver.iterations);
fprintf('SVDs       : %d\n', out_fbs.solver.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.solver.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.solver.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('\n');
fprintf('iterations : %d\n', out_lbfgs.solver.iterations);
fprintf('SVDs       : %d\n', out_lbfgs.solver.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.solver.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.solver.residual(end));
