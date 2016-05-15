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

x_star = [0; 0; 2.114635341704963e-01; 0; 2.845881348733116e+00];

ASSERT_TOL = 1e-8;

%% adaptive

baseopt.display = 0;
baseopt.adaptive = 1;
baseopt.tol = 1e-14;
baseopt.maxit = 10000;

opt_fbs = baseopt; opt_fbs.solver = 'fbs'; opt_fbs.variant = 'basic';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);

assert(out_fbs.iterations < baseopt.maxit);

opts = {};
outs = {};

opts{end+1} = baseopt; opts{end}.solver = 'fbs'; opts{end}.variant = 'fast'; 
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, x0, aff, [], opts{i});
    assert(outs{i}.iterations < opts{i}.maxit);
    assert(outs{i}.objective(end) - out_fbs.objective(end) <= ASSERT_TOL);
    assert(norm(outs{i}.x - x_star, 'fro') <= ASSERT_TOL);
end
