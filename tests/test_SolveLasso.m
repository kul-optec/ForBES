A = [1,  2, -1, -1; ...
    -2, -1,  0, -1; ...
    3,  0,  4, -1; ...
    -4, -1, -3,  1; ...
    5,  3,  2,  3]';
b = [1, 2, 3, 4]';

[m, n] = size(A);

f = quadLoss(1, zeros(m,1));
aff = {A, -b};
lam = 5.0;
g = l1Norm(lam);
x0 = zeros(n, 1);

%% non-adaptive

baseopt.display = 0;
baseopt.adaptive = 0;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic'; opt_fbs.tol = 1e-12;
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);

opt_afbs = baseopt; opt_afbs.solver = 'fbs'; opt_afbs.variant = 'fast'; opt_afbs.tol = 1e-12;
out_afbs = forbes(f, g, x0, aff, [], opt_afbs);

assert(norm(out_fbs.x - out_afbs.x) <= 1e-11);

opt_minfbe_bfgs = baseopt; opt_minfbe_bfgs.solver = 'minfbe'; opt_minfbe_bfgs.method = 'bfgs'; opt_minfbe_bfgs.tol = 1e-12;
out_minfbe_bfgs = forbes(f, g, x0, aff, [], opt_minfbe_bfgs);

assert(norm(out_fbs.x - out_minfbe_bfgs.x) <= 1e-11);

opt_minfbe_lbfgs = baseopt; opt_minfbe_lbfgs.solver = 'minfbe'; opt_minfbe_lbfgs.method = 'lbfgs'; opt_minfbe_lbfgs.tol = 1e-12;
out_minfbe_lbfgs = forbes(f, g, x0, aff, [], opt_minfbe_lbfgs);

assert(norm(out_fbs.x - out_minfbe_lbfgs.x) <= 1e-11);

opt_zerofpr_bfgs = baseopt; opt_zerofpr_bfgs.solver = 'zerofpr'; opt_zerofpr_bfgs.method = 'bfgs'; opt_zerofpr_bfgs.tol = 1e-12;
out_zerofpr_bfgs = forbes(f, g, x0, aff, [], opt_zerofpr_bfgs);

assert(norm(out_fbs.x - out_zerofpr_bfgs.x) <= 1e-11);

opt_zerofpr_lbfgs = baseopt; opt_zerofpr_lbfgs.solver = 'zerofpr'; opt_zerofpr_lbfgs.method = 'lbfgs'; opt_zerofpr_lbfgs.tol = 1e-12;
out_zerofpr_lbfgs = forbes(f, g, x0, aff, [], opt_zerofpr_lbfgs);

assert(norm(out_fbs.x - out_zerofpr_lbfgs.x) <= 1e-11);

%% non-adaptive

baseopt.display = 0;
baseopt.adaptive = 1;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic'; opt_fbs.tol = 1e-12;
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);

opt_afbs = baseopt; opt_afbs.solver = 'fbs'; opt_afbs.variant = 'fast'; opt_afbs.tol = 1e-12;
out_afbs = forbes(f, g, x0, aff, [], opt_afbs);

assert(norm(out_fbs.x - out_afbs.x) <= 1e-11);

opt_minfbe_bfgs = baseopt; opt_minfbe_bfgs.solver = 'minfbe'; opt_minfbe_bfgs.method = 'bfgs'; opt_minfbe_bfgs.tol = 1e-12;
out_minfbe_bfgs = forbes(f, g, x0, aff, [], opt_minfbe_bfgs);

assert(norm(out_fbs.x - out_minfbe_bfgs.x) <= 1e-11);

opt_minfbe_lbfgs = baseopt; opt_minfbe_lbfgs.solver = 'minfbe'; opt_minfbe_lbfgs.method = 'lbfgs'; opt_minfbe_lbfgs.tol = 1e-12;
out_minfbe_lbfgs = forbes(f, g, x0, aff, [], opt_minfbe_lbfgs);

assert(norm(out_fbs.x - out_minfbe_lbfgs.x) <= 1e-11);

opt_zerofpr_bfgs = baseopt; opt_zerofpr_bfgs.solver = 'zerofpr'; opt_zerofpr_bfgs.method = 'bfgs'; opt_zerofpr_bfgs.tol = 1e-12;
out_zerofpr_bfgs = forbes(f, g, x0, aff, [], opt_zerofpr_bfgs);

assert(norm(out_fbs.x - out_zerofpr_bfgs.x) <= 1e-11);

opt_zerofpr_lbfgs = baseopt; opt_zerofpr_lbfgs.solver = 'zerofpr'; opt_zerofpr_lbfgs.method = 'lbfgs'; opt_zerofpr_lbfgs.tol = 1e-12;
out_zerofpr_lbfgs = forbes(f, g, x0, aff, [], opt_zerofpr_lbfgs);

assert(norm(out_fbs.x - out_zerofpr_lbfgs.x) <= 1e-11);