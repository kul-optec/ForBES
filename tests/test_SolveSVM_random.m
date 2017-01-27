close all;
clear;

% solve a small SVM problem via the dual QP using quadprog
% then compare with the solution found with ForBES

n = 2100; % number of features (= number of variables minus one)
m = 130; % number of samples

w = randn(n, 1);  % N(0,1), 30% sparse

A = randn(m, n);
btrue = sign(A*w);

% noise is function of problem size use 0.1 for large problem
b = sign(btrue + sqrt(0.1)*randn(m,1)); % labels with noise
mu = 1.0;

% solve dual problem using QUADPROG
BA = diag(sparse(b))*A;
Q = BA*BA';
q = ones(size(A, 1), 1);
opt_qp = optimoptions('quadprog','Display','off');
[lambda_qp, fval_qp, flag_qp, output_qp] = quadprog(Q, -q, [], [], [], [], 0, mu, [], opt_qp);
x_qp = BA'*lambda_qp;

f = quadLoss();
g = hingeLoss(mu, b);
constr = {A, -1, zeros(m, 1)};
y0 = zeros(m, 1);

TOL = 1e-12;
ASSERT_TOLX = 1e-6;
ASSERT_TOLF = 1e-8;

%% adaptive

baseopt.display = 0;
baseopt.adaptive = 1;
baseopt.maxit = 10000;
baseopt.tol = TOL;

opts = {};
outs = {};

opts{end+1} = baseopt; opts{end}.solver = 'fbs'; opts{end}.variant = 'basic';
opts{end+1} = baseopt; opts{end}.solver = 'fbs'; opts{end}.variant = 'fast';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking-armijo';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, y0, [], constr, opts{i});
    assert(outs{i}.flag == 0);
    assert(norm(A*outs{i}.x1 - outs{i}.z, 'inf') <= 10*TOL);
    assert(abs(outs{i}.dual.objective(end) - fval_qp)/(1+abs(fval_qp)) <= ASSERT_TOLF);
    assert(norm(outs{i}.x1 - x_qp, inf)/(1+norm(x_qp, inf)) <= ASSERT_TOLX);
    fprintf('.');
end
