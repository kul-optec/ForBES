clear;

A = [1,  2, -1, -1; ...
    -2, -1,  0, -1; ...
    3,  0,  4, -1; ...
    -4, -1, -3,  1; ...
    5,  3,  2,  3]';
b = [1, 2, 3, 4]';

[m, n] = size(A);

f = logLoss(1.0);
aff = {A, -b};
lam = 0.1;
g = l1Norm(lam);
x0 = zeros(n, 1);

ASSERT_TOL = 1e-6;

%% adaptive

baseopt.display = 0;
baseopt.adaptive = 1;
baseopt.tol = 1e-14;
baseopt.maxit = 10000;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);

assert(out_fbs.iterations < baseopt.maxit);

opt_afbs = baseopt; opt_afbs.solver = 'fbs'; opt_afbs.variant = 'fast';
out_afbs = forbes(f, g, x0, aff, [], opt_afbs);

assert(norm(out_fbs.x - out_afbs.x, 'inf') <= ASSERT_TOL);

opt_minfbe_g_bfgs = baseopt; opt_minfbe_g_bfgs.solver = 'minfbe'; opt_minfbe_g_bfgs.method = 'bfgs';
out_minfbe_g_bfgs = forbes(f, g, x0, aff, [], opt_minfbe_g_bfgs);

assert(norm(out_fbs.x - out_minfbe_g_bfgs.x, 'inf') <= ASSERT_TOL);

opt_minfbe_g_lbfgs = baseopt; opt_minfbe_g_lbfgs.solver = 'minfbe'; opt_minfbe_g_lbfgs.method = 'lbfgs';
out_minfbe_g_lbfgs = forbes(f, g, x0, aff, [], opt_minfbe_g_lbfgs);

assert(norm(out_fbs.x - out_minfbe_g_lbfgs.x, 'inf') <= ASSERT_TOL);

opt_zerofpr_bfgs = baseopt; opt_zerofpr_bfgs.solver = 'zerofpr'; opt_zerofpr_bfgs.method = 'bfgs';
out_zerofpr_bfgs = forbes(f, g, x0, aff, [], opt_zerofpr_bfgs);

assert(norm(out_fbs.x - out_zerofpr_bfgs.x, 'inf') <= ASSERT_TOL);

opt_zerofpr_lbfgs = baseopt; opt_zerofpr_lbfgs.solver = 'zerofpr'; opt_zerofpr_lbfgs.method = 'lbfgs';
out_zerofpr_lbfgs = forbes(f, g, x0, aff, [], opt_zerofpr_lbfgs);

assert(norm(out_fbs.x - out_zerofpr_lbfgs.x, 'inf') <= ASSERT_TOL);

opt_minfbe_b_bfgs = baseopt; opt_minfbe_b_bfgs.solver = 'minfbe'; opt_minfbe_b_bfgs.variant = 'basic'; opt_minfbe_b_bfgs.method = 'bfgs';
out_minfbe_b_bfgs = forbes(f, g, x0, aff, [], opt_minfbe_b_bfgs);

assert(norm(out_fbs.x - out_minfbe_b_bfgs.x, 'inf') <= ASSERT_TOL);

opt_minfbe_b_lbfgs = baseopt; opt_minfbe_b_lbfgs.solver = 'minfbe'; opt_minfbe_b_lbfgs.variant = 'basic'; opt_minfbe_b_lbfgs.method = 'lbfgs';
out_minfbe_b_lbfgs = forbes(f, g, x0, aff, [], opt_minfbe_b_lbfgs);

assert(norm(out_fbs.x - out_minfbe_b_lbfgs.x, 'inf') <= ASSERT_TOL);