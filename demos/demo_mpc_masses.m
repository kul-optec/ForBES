close all;
clear;

rng(0);

S = system_masses(3);

% Get dynamics from the system
A = S.A;
B = S.B;
[n, m] = size(B);

% Define cost matrices
Q = 1*eye(n);
R = 1*eye(m);
[~, Q_f] = dlqr(A,B,Q,R);

% Input constraints
umin = -10*ones(m,1);
umax = +10*ones(m,1);
xmin = -1*ones(n,1);
xmax = +1*ones(n,1);
% Initial point and prediction horizon
x0 = randn(n, 1);
N = 20;

% Define ForBES problem
f = lqrCost(x0, Q, R, Q_f, A, B, N);

% one block constraints
lb = [repmat([xmin;umin],N,1);xmin];
ub = [repmat([xmax;umax],N,1);xmax];
weights = [repmat([1e2*ones(n,1);inf*ones(m,1)],N,1);1e2*ones(n,1)];
g = distBox(lb,ub,weights);

% % many blocks constraints
% prob.g = separableSum(  {distBox(xmin,xmax,1e5*ones(n,1)),...
%     indBox(umin,umax), distBox(-ones(n,1),ones(n,1),1e5*ones(n,1))}, ...
%     [repmat([1, 2], 1, N), 3], ...    % indices
%     [repmat([n, m], 1, N), n]);       % dimensions

% % two blocks constraints
% lb = repmat([xmin;umin],N,1);
% ub = repmat([xmax;umax],N,1);
% weights = repmat([1e5*ones(n,1);inf*ones(m,1)],N,1);
% prob.g = separableSum(  {distBox(lb,ub,weights), distBox(-ones(n,1),ones(n,1),1e5*ones(n,1))}, ...
%     [1, 2], ...    % indices
%     [N*(n+m), n]);       % dimensions

b = zeros(n*(N+1)+m*N,1);
constr = {1, -1, b};
y0 = zeros(n*(N+1)+m*N,1);
opt.maxit = 10000;
opt.tol = 1e-9;
opt.display = 1;

fprintf('\nFast FBS\n');
opt_fbs = opt;
opt_fbs.method = 'fbs';
opt_fbs.variant = 'fast';
out = forbes(f, g, y0, [], constr, opt_fbs);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('matvecs    : %d\n', out.operations.C1);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nL-BFGS\n');
opt_lbfgs = opt;
opt_lbfgs.method = 'lbfgs';
out = forbes(f, g, y0, [], constr, opt_lbfgs);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('matvecs    : %d\n', out.operations.C1);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

fprintf('\nCG-DYHS\n');
opt_cg = opt;
opt_cg.method = 'cg-dyhs';
out = forbes(f, g, y0, [], constr, opt_cg);
fprintf('message    : %s\n', out.message);
fprintf('iterations : %d\n', out.iterations);
fprintf('matvecs    : %d\n', out.operations.C1);
fprintf('time       : %7.4e\n', out.ts(end));
fprintf('residual   : %7.4e\n', out.residual(end));

