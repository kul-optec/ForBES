% solve a sparse logistic regression problem using ForBES

close all;
clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 6000;
n = 50000;
x_orig = sprandn(n, 1, 30/n);
A = sprandn(m, n, 200/n);
b = 2*(rand(m,1) <= 1./(1+exp(-A*x_orig))) - 1;
lam_max = norm(0.5*(A'*b),'inf')/m;
lam = 0.5*lam_max;

f = logLoss(1/m);
aff = {diag(sparse(b))*A};
g = l1Norm(lam);
x0 = zeros(n, 1);
opt.maxit = 1000;
opt.tol = 1e-12;
opt.adaptive = 1;
opt.display = 2;

fprintf('\nFast FBS\n');
opt_fbs = opt;
opt_fbs.method = 'fbs';
opt_fbs.variant = 'fast';
out = forbes(f, g, x0, aff, [], opt_fbs);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('f evals    : %d\n', out.operations.cnt_f2);
fprintf('matvecs    : %d\n', out.operations.cnt_C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out = forbes(f, g, x0, aff, [], opt_lbfgs);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('f evals    : %d\n', out.operations.cnt_f2);
fprintf('matvecs    : %d\n', out.operations.cnt_C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nCG-DYHS\n');
opt_cg1 = opt;
opt_cg1.method = 'cg-dyhs';
out = forbes(f, g, x0, aff, [], opt_cg1);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('f evals    : %d\n', out.operations.cnt_f2);
fprintf('matvecs    : %d\n', out.operations.cnt_C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nCG-PRP\n');
opt_cg2 = opt;
opt_cg2.method = 'cg-prp';
out = forbes(f, g, x0, aff, [], opt_cg2);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('f evals    : %d\n', out.operations.cnt_f2);
fprintf('matvecs    : %d\n', out.operations.cnt_C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nCG-DESCENT\n');
opt_cg3 = opt;
opt_cg3.method = 'cg-desc';
out = forbes(f, g, x0, aff, [], opt_cg3);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('f evals    : %d\n', out.operations.cnt_f2);
fprintf('matvecs    : %d\n', out.operations.cnt_C2);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));
