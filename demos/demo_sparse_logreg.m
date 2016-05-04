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
aff = {diag(sparse(b))*A};
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
fprintf('iterations : %d\n', out_fbs.iterations);
fprintf('f evals    : %d\n', out_fbs.operations.f2);
fprintf('f'' evals   : %d\n', out_fbs.operations.gradf2);
fprintf('matvecs    : %d\n', out_fbs.operations.C2);
fprintf('g          : %d\n', out_fbs.operations.g);
fprintf('prox       : %d\n', out_fbs.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('message    : %s\n', out_lbfgs.message);
fprintf('iterations : %d\n', out_lbfgs.iterations);
fprintf('f evals    : %d\n', out_lbfgs.operations.f2);
fprintf('f'' evals   : %d\n', out_lbfgs.operations.gradf2);
fprintf('matvecs    : %d\n', out_lbfgs.operations.C2);
fprintf('g          : %d\n', out_lbfgs.operations.g);
fprintf('prox       : %d\n', out_lbfgs.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.residual(end));

fprintf('\nCG-DYHS\n');
opt_cg1 = opt;
opt_cg1.method = 'cg-dyhs';
out_cg1 = forbes(f, g, x0, aff, [], opt_cg1);
fprintf('message    : %s\n', out_cg1.message);
fprintf('iterations : %d\n', out_cg1.iterations);
fprintf('f evals    : %d\n', out_cg1.operations.f2);
fprintf('f'' evals   : %d\n', out_cg1.operations.gradf2);
fprintf('matvecs    : %d\n', out_cg1.operations.C2);
fprintf('g          : %d\n', out_cg1.operations.g);
fprintf('prox       : %d\n', out_cg1.operations.proxg);
fprintf('time       : %7.4e\n', out_cg1.ts(end));
fprintf('residual   : %7.4e\n', out_cg1.residual(end));

