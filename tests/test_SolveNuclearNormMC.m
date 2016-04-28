clear;

rng(0, 'twister'); % uncomment this to control the random number generator

m = 50; % number of rows
n = 50; % number of column of the original matrix M
d = 1000/(m*n); % density of coefficients sampled from M
r = 3; % rank of M

U = randn(m, r);
V = randn(n, r);
M = U*V';

P = sprand(m, n, d) ~= 0; % sampling pattern
B = full(M.*P);

f = quadLoss(P(:), B(:));
lam = 5;
g = nuclearNorm(m, n, lam, 'inexact');
x0 = zeros(m*n, 1);

ASSERT_TOL = 1e-8;

%% run methods

baseopt.display = 0;
baseopt.adaptive = 0;
baseopt.tol = 1e-10;
baseopt.maxit = 1000;
baseopt.Lf = 1;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, x0, [], [], opt_fbs);

assert(out_fbs.iterations < baseopt.maxit);

opt_afbs = baseopt; opt_afbs.solver = 'fbs'; opt_afbs.variant = 'fast';
out_afbs = forbes(f, g, x0, [], [], opt_afbs);

assert(norm(out_fbs.x - out_afbs.x) <= ASSERT_TOL);

opt_minfbe_g_lbfgs = baseopt; opt_minfbe_g_lbfgs.solver = 'minfbe'; opt_minfbe_g_lbfgs.method = 'lbfgs';
out_minfbe_g_lbfgs = forbes(f, g, x0, [], [], opt_minfbe_g_lbfgs);

assert(norm(out_fbs.x - out_minfbe_g_lbfgs.x) <= ASSERT_TOL);

opt_zerofpr_lbfgs = baseopt; opt_zerofpr_lbfgs.solver = 'zerofpr'; opt_zerofpr_lbfgs.method = 'lbfgs';
out_zerofpr_lbfgs = forbes(f, g, x0, [], [], opt_zerofpr_lbfgs);

assert(norm(out_fbs.x - out_zerofpr_lbfgs.x) <= ASSERT_TOL);

opt_minfbe_b_lbfgs = baseopt; opt_minfbe_b_lbfgs.solver = 'minfbe'; opt_minfbe_b_lbfgs.variant = 'basic'; opt_minfbe_b_lbfgs.method = 'lbfgs';
out_minfbe_b_lbfgs = forbes(f, g, x0, [], [], opt_minfbe_b_lbfgs);

assert(norm(out_fbs.x - out_minfbe_b_lbfgs.x) <= ASSERT_TOL);
