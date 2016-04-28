clear;

rng(0, 'twister');

n = 2000; % number of features (= number of variables minus one)
m = 300; % number of samples

w = sprandn(n, 1, 0.3);  % N(0,1), 30% sparse
v = randn(1);            % random intercept

X = sprandn(m, n, 10/n);
btrue = sign(X*w + v);

% noise is function of problem size use 0.1 for large problem
b = sign(X*w + v + sqrt(0.1)*randn(m,1)); % labels with noise

A = [X, ones(m, 1)];

ratio = sum(b == 1)/(m);
lam = 0.1 * norm((1-ratio)*sum(A(b==1,:),1) + ratio*sum(A(b==-1,:),1), 'inf');

f = quadLoss(lam, zeros(n+1, 1));
g = hingeLoss(1, b);
constr = {A, -1, zeros(m, 1)};
y0 = zeros(m, 1);

ASSERT_TOL = 1e-10;

%% adaptive

baseopt.display = 0;
baseopt.adaptive = 1;
baseopt.maxit = 10000;
baseopt.tol = 1e-14;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, y0, [], constr, opt_fbs);

assert(out_fbs.iterations < baseopt.maxit);

opt_afbs = baseopt; opt_afbs.solver = 'fbs'; opt_afbs.variant = 'fast';
out_afbs = forbes(f, g, y0, [], constr, opt_afbs);

assert(norm(out_fbs.y - out_afbs.y, 'inf') <= ASSERT_TOL);

opt_minfbe_g_bfgs = baseopt; opt_minfbe_g_bfgs.solver = 'minfbe'; opt_minfbe_g_bfgs.method = 'bfgs';
out_minfbe_g_bfgs = forbes(f, g, y0, [], constr, opt_minfbe_g_bfgs);

assert(norm(out_fbs.y - out_minfbe_g_bfgs.y, 'inf') <= ASSERT_TOL);

opt_minfbe_g_lbfgs = baseopt; opt_minfbe_g_lbfgs.solver = 'minfbe'; opt_minfbe_g_lbfgs.method = 'lbfgs';
out_minfbe_g_lbfgs = forbes(f, g, y0, [], constr, opt_minfbe_g_lbfgs);

assert(norm(out_fbs.y - out_minfbe_g_lbfgs.y, 'inf') <= ASSERT_TOL);

opt_zerofpr_bfgs = baseopt; opt_zerofpr_bfgs.solver = 'zerofpr'; opt_zerofpr_bfgs.method = 'bfgs';
out_zerofpr_bfgs = forbes(f, g, y0, [], constr, opt_zerofpr_bfgs);

assert(norm(out_fbs.y - out_zerofpr_bfgs.y, 'inf') <= ASSERT_TOL);

opt_zerofpr_lbfgs = baseopt; opt_zerofpr_lbfgs.solver = 'zerofpr'; opt_zerofpr_lbfgs.method = 'lbfgs';
out_zerofpr_lbfgs = forbes(f, g, y0, [], constr, opt_zerofpr_lbfgs);

assert(norm(out_fbs.y - out_zerofpr_lbfgs.y, 'inf') <= ASSERT_TOL);

opt_minfbe_b_bfgs = baseopt; opt_minfbe_b_bfgs.solver = 'minfbe'; opt_minfbe_b_bfgs.variant = 'basic'; opt_minfbe_b_bfgs.method = 'bfgs';
out_minfbe_b_bfgs = forbes(f, g, y0, [], constr, opt_minfbe_b_bfgs);

assert(norm(out_fbs.y - out_minfbe_b_bfgs.y, 'inf') <= ASSERT_TOL);

opt_minfbe_b_lbfgs = baseopt; opt_minfbe_b_lbfgs.solver = 'minfbe'; opt_minfbe_b_lbfgs.variant = 'basic'; opt_minfbe_b_lbfgs.method = 'lbfgs';
out_minfbe_b_lbfgs = forbes(f, g, y0, [], constr, opt_minfbe_b_lbfgs);

assert(norm(out_fbs.y - out_minfbe_b_lbfgs.y, 'inf') <= ASSERT_TOL);