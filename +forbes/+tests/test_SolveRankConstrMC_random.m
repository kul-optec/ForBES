close all;
clear;

% rng(0, 'twister'); % uncomment this to control the random number generator

m = 30; % number of rows
n = 30; % number of column of the original matrix M
d = 0.5; % density of coefficients sampled from M
r = 5; % rank of M
r_target = 3;

U = randn(m, r);
V = randn(n, r);
M = U*V';

P = sprand(m, n, d) ~= 0; % sampling pattern
B = full(M.*P);

f = forbes.functions.SqrNormL2(P);
aff = {1, -B};
g = forbes.functions.IndRankBall(r_target);
x0 = zeros(m, n);

ASSERT_TOL = 1e-5;

%% run methods

baseopt.display = 0;
baseopt.tol = 1e-6;
baseopt.maxit = 1000;
baseopt.Lf = 1;
baseopt.report = 1;

opts = {};
outs = {};

opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking-nm';
% opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking';
% opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking-nm';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking-nm';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, x0, aff, [], opts{i});
    assert(outs{i}.solver.iterations < opts{i}.maxit);
    assert(norm(outs{i}.solver.residual(end), 'inf') <= ASSERT_TOL);
    fprintf('.');
end
