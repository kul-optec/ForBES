% solve a nuclear norm regularized matrix completion problem using ForBES

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

f = forbes.functions.SqrNormL2(P);
aff = {1, -B};
g = forbes.functions.NuclearNorm(lam, 'inexact');
x0 = zeros(m, n);
opt.maxit = 1000;
opt.tol = 1e-6;
opt.Lf = 1;
opt.display = 1;

fprintf('\nFast FBS\n');
opt_fbs = opt;
opt_fbs.solver = 'fbs';
opt_fbs.variant = 'fast';
out = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('\n');
fprintf('iterations : %d\n', out.solver.iterations);
fprintf('SVDs       : %d\n', out.solver.operations.proxg);
fprintf('time       : %7.4e\n', out.solver.ts(end));
fprintf('residual   : %7.4e\n', out.solver.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('\n');
fprintf('iterations : %d\n', out.solver.iterations);
fprintf('SVDs       : %d\n', out.solver.operations.proxg);
fprintf('time       : %7.4e\n', out.solver.ts(end));
fprintf('residual   : %7.4e\n', out.solver.residual(end));
