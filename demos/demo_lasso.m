% solve a lasso problem using ForBES

close all;
clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 5000; % number of observations
n = 30000; % number of features
x_orig = sprandn(n, 1, 20/n); % generate random sparse model
A = sprandn(m, n, 100/n); % generate random sparse design matrix
b = A*x_orig + randn(m, 1)/10; % compute labels and add noise

fprintf('%d nonzero features\n', nnz(A));
fprintf('%.2f nnz per row\n', nnz(A)/numel(A)*n);

% for lam >= lam_max the solution is zero
lam_max = norm(A'*b,'inf');
lam = 0.05*lam_max;

f = forbes.functions.SqrNormL2();
aff = {A, -b};
g = forbes.functions.NormL1(lam);
x0 = zeros(n, 1);
opt.maxit = 10000;
opt.tol = 1e-9;
opt.display = 1;

fprintf('\nFast FBS\n');
opt_fbs = opt;
opt_fbs.solver = 'fbs';
opt_fbs.variant = 'fast';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('\n');
fprintf('message    : %s\n', out_fbs.message);
fprintf('iterations : %d\n', out_fbs.solver.iterations);
fprintf('matvecs    : %d\n', out_fbs.solver.operations.C1);
fprintf('prox       : %d\n', out_fbs.solver.operations.proxg);
fprintf('time       : %7.4e\n', out_fbs.solver.ts(end));
fprintf('residual   : %7.4e\n', out_fbs.solver.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out_lbfgs = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('\n');
fprintf('message    : %s\n', out_lbfgs.message);
fprintf('iterations : %d\n', out_lbfgs.solver.iterations);
fprintf('matvecs    : %d\n', out_lbfgs.solver.operations.C1);
fprintf('prox       : %d\n', out_lbfgs.solver.operations.proxg);
fprintf('time       : %7.4e\n', out_lbfgs.solver.ts(end));
fprintf('residual   : %7.4e\n', out_lbfgs.solver.residual(end));
