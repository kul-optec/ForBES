close all;
clear;

A = [1,  2, -1, -1; ...
    -2, -1,  0, -1; ...
    3,  0,  4, -1; ...
    -4, -1, -3,  1; ...
    5,  3,  2,  3]';
b = [1, 2, 3, 4]';

[m, n] = size(A);

f = forbes.functions.SqrNormL2();
aff = {A, -b};
lam = 0.1*norm(A'*b, 'inf');
g = forbes.functions.NormL1(lam);
x0 = zeros(n, 1);

x_star = [-3.877278911564627e-01; 0; 0; 2.174149659863943e-02; 6.168435374149660e-01];

ASSERT_TOL = 1e-6;

baseopt.display = 2;
baseopt.tol = 1e-8;
baseopt.maxit = 10000;

opt_fbs = baseopt; opt_fbs.solver = 'fbs';
out_fbs = forbes(f, g, x0, aff, [], opt_fbs);

assert(out_fbs.solver.iterations < baseopt.maxit);

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
opts{end+1} = baseopt; opts{end}.solver = 'nama'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'nama'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'nama'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'nama'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'nama'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, x0, aff, [], opts{i});
    assert(outs{i}.flag == 0);
    assert(norm(outs{i}.x - out_fbs.x,inf)/(1+norm(out_fbs.x,inf)) <= ASSERT_TOL);
    fprintf('.');
end
