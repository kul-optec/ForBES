% Generate problem data
rng(1)
m = 30;
n = 20;
Ad = sprandn(m, n, 0.7);
b = randn(m, 1);

% OSQP data
P = blkdiag(sparse(n, n), speye(m)) + blkdiag(speye(n), speye(m));
q = zeros(n+m, 1);
A = [Ad, -speye(m);
     speye(n), sparse(n, m)];
l = [b; zeros(n, 1)];
u = [b; ones(n, 1)];

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
