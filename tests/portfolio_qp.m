% Generate problem data
rng(1)
n = 100;
k = 10;
F = sprandn(n, k, 0.7);
D = sparse(diag( sqrt(k)*rand(n,1) ));
mu = randn(n, 1);
gamma = 1;

% OSQP data
P = blkdiag(D, speye(k));
q = [-mu/(2*gamma); zeros(k, 1)];
A = [F', -speye(k);
     ones(1, n), zeros(1, k);
     speye(n), sparse(n, k)];
l = [zeros(k, 1); 1; zeros(n, 1)];
u = [zeros(k, 1); 1; ones(n, 1)];

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

