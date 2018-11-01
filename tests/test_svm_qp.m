% Generate problem data
rng(1)
n = 10;
m = 1000;
N = ceil(m/2);
gamma = 1;
A_upp = sprandn(N, n, 0.5);
A_low = sprandn(N, n, 0.5);
Ad = [A_upp / sqrt(n) + (A_upp ~= 0) / n;
      A_low / sqrt(n) - (A_low ~= 0) / n];
b = [ones(N, 1); -ones(N,1)];

% OSQP data
P = blkdiag(speye(n), sparse(m, m)) + blkdiag(speye(n), speye(m));
q = [zeros(n,1); gamma*ones(m,1)];
A = [diag(b)*Ad, -speye(m);
     sparse(m, n), speye(m)];
l = [-inf*ones(m, 1); zeros(m, 1)];
u = [-ones(m, 1); inf*ones(m, 1)];

%% QUADPROG SANITY CHECK

Aqp = [A;-A];
bqp = [u;-l];
tic
[x,f,flag] = quadprog(P,q,Aqp,bqp);
elapsedTime1 = toc;



%% Create an OSQP object
prob = osqp;

% Setup workspace
prob.setup(P, q, A, l, u);

% Solve problem
tic,
out = prob.solve();
elapsedTime2 = toc;

%% NAMA Solver Test
opt.solver = 'nama';
opt.term = @(prob, it, cache_0, cache_x,scale) custom_termination(prob, it, cache_0, cache_x,scale);
tic,
nama_out = forbes_qp(P, q, A, l, u,[],[], [], [], opt);
elapsedTime3 = toc;

fprintf('\nThe time for the solvers\n');
fprintf('quadprog: %.2e, osqp:%.2e, nama:%.2e\n', elapsedTime1, elapsedTime2, elapsedTime3);
fprintf('osqp:(%.2e,%2d), nama:(%.2e,%2d)\n', out.info.pri_res, out.info.iter, nama_out.solver.solver.residual(end),nama_out.solver.solver.iterations);

return;

