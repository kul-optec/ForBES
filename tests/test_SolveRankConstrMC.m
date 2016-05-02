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

ASSERT_TOL = 1e-10;

%% run methods

baseopt.display = 0;
baseopt.adaptive = 0;
baseopt.tol = 1e-10;
baseopt.maxit = 1000;
baseopt.Lf = 1;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, x0, [], [], opt_fbs);

assert(out_fbs.iterations < baseopt.maxit);

opts = {};
outs = {};

opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking-nm';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, x0, [], [], opts{i});
    assert(outs{i}.iterations < opts{i}.maxit);
    assert(norm(outs{i}.residual(end), 'inf') <= ASSERT_TOL);
end
