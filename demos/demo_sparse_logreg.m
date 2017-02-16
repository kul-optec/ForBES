% solve a sparse logistic regression problem using ForBES

close all;
clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 2000;
n = 50000;
x_orig = sprandn(n, 1, 20/n);
A = sprandn(m, n, 500/n);
b = 2*(rand(m,1) <= 1./(1+exp(-A*x_orig))) - 1;
lam_max = norm(0.5*(A'*b),'inf')/m;
lam = 0.1*lam_max;

f = logLoss(1/m);
aff = {diag(sparse(b))*A, zeros(m, 1)};
g = l1Norm(lam);
x0 = zeros(n, 1);
opt.maxit = 10000;
opt.tol = 1e-6;
opt.adaptive = 1;
opt.display = 1;

fprintf('\nFast FBS\n');
opt_fbs = opt;
opt_fbs.solver = 'fbs';
opt_fbs.variant = 'fast';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('message    : %s\n', out_fbs.message);
fprintf('iterations : %d\n', out_fbs.solver.iterations);
fprintf('f evals    : %d\n', out_fbs.solver.operations.f2);
fprintf('f'' evals   : %d\n', out_fbs.solver.operations.gradf2);
fprintf('matvecs    : %d\n', out_fbs.solver.operations.C2);
fprintf('g          : %d\n', out_fbs.solver.operations.g);
fprintf('prox       : %d\n', out_fbs.solver.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.solver.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.solver.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('message    : %s\n', out_lbfgs.message);
fprintf('iterations : %d\n', out_lbfgs.solver.iterations);
fprintf('f evals    : %d\n', out_lbfgs.solver.operations.f2);
fprintf('f'' evals   : %d\n', out_lbfgs.solver.operations.gradf2);
fprintf('matvecs    : %d\n', out_lbfgs.solver.operations.C2);
fprintf('g          : %d\n', out_lbfgs.solver.operations.g);
fprintf('prox       : %d\n', out_lbfgs.solver.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.solver.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.solver.residual(end));
