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
g = indRankBall(m, n, r);
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

opt_zerofpr_lbfgs = baseopt; opt_zerofpr_lbfgs.solver = 'zerofpr'; opt_zerofpr_lbfgs.method = 'lbfgs';
out_zerofpr_lbfgs = forbes(f, g, x0, [], [], opt_zerofpr_lbfgs);

assert(norm(out_fbs.x - out_zerofpr_lbfgs.x, 'inf') <= ASSERT_TOL);
