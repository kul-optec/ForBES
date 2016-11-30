% rng(0, 'twister');

m = 50;
n = 200;

A = randn(m, n);
b = randn(m, 1);

f = quadLoss();
aff = {A, -b};
lam = 0.3*norm(A'*b, 'inf');
g = l1Norm(lam);
x0 = zeros(n, 1);

ASSERT_TOL = 1e-8;

baseopt.display = 0;
baseopt.tol = 1e-12;
baseopt.maxit = 10000;

opt_fbs = baseopt; opt_fbs.solver = 'fbs';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);

assert(out_fbs.iterations < baseopt.maxit);

opts = {};
outs = {};

opts{end+1} = baseopt; opts{end}.solver = 'fbs'; opts{end}.variant = 'fast';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'broyden';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, x0, aff, [], opts{i});
    assert(outs{i}.iterations < opts{i}.maxit);
    assert(norm(outs{i}.x - out_fbs.x,inf)/(1+norm(out_fbs.x,inf)) <= ASSERT_TOL);
    fprintf('.');
end
