close all;
clear;

n = 20;
m = 60;
densQ = 0.2;
densA = 0.2;
act = 3;

Q = sprandsym(n, densQ, 1, 1)+1e-1*speye(n);
A = sprandn(m, n, densA);
x_star = randn(n, 1);
y_star = [rand(act, 1); zeros(m-act, 1)];
q = -Q*x_star - A'*y_star;
b = [A(1:act,:)*x_star; A(act+1:end,:)*x_star + rand(m-act,1)];
f_star = 0.5*(x_star'*(Q*x_star)) + q'*x_star;

f = quadratic(Q, q);
g = indPos();
constr = {A, 1, b};
y0 = zeros(m, 1);

TOL = 1e-8;
ASSERT_TOLX = 1e-6;
ASSERT_TOLF = 1e-10;

% run solvers

baseopt.display = 0;
baseopt.adaptive = 0;
baseopt.maxit = 10000;
baseopt.tol = TOL;

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
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'bfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'broyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'lbroyden'; opts{end}.linesearch = 'backtracking';
opts{end+1} = baseopt; opts{end}.solver = 'amls'; opts{end}.method = 'rbroyden'; opts{end}.linesearch = 'backtracking';

for i = 1:length(opts)
    outs{end+1} = forbes(f, g, y0, [], constr, opts{i});
    assert(outs{i}.iterations < opts{i}.maxit);
    assert(norm(A*outs{i}.x1 + outs{i}.z - b, 'inf') <= 10*TOL);
    assert(abs(outs{i}.dual.objective(end) + f_star)/(1+abs(f_star)) <= ASSERT_TOLF);
    assert(norm(outs{i}.x1 - x_star, inf)/(1+norm(x_star, inf)) <= ASSERT_TOLX);
    fprintf('.');
end
