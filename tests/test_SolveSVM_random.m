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

opts = {};
outs = {};

opts{end+1} = baseopt; opts{end}.solver = 'fbs'; opts{end}.variant = 'fast';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, y0, [], constr, opts{i});
    assert(outs{i}.iterations < opts{i}.maxit);
    assert(norm(outs{i}.y - out_fbs.y, 'inf') <= ASSERT_TOL);
end
