% Generate problem data
rng(1)
n = 10;
m = 1000;
Ad = sprandn(m, n, 0.5);
x_true = (randn(n, 1) > 0.8) .* randn(n, 1) / sqrt(n);
b = Ad * x_true + 0.5 * randn(m, 1);
gammas = linspace(1, 10, 11);

% OSQP data
P = blkdiag(sparse(n, n), speye(m), sparse(n, n)) + blkdiag(speye(n), speye(m), speye(n));
q = zeros(2*n+m, 1);
A = [Ad, -speye(m), sparse(m,n);
    speye(n), sparse(n, m), -speye(n);
    speye(n), sparse(n, m), speye(n);];
l = [b; -inf*ones(n, 1); zeros(n, 1)];
u = [b; zeros(n, 1); inf*ones(n, 1)];

%% QUADPROG SANITY CHECK

Aqp = [A;-A];
bqp = [u;-l];
tic
[x,f,flag] = quadprog(P,q,Aqp,bqp);
elapsedTime1 = toc;


%% Create an OSQP object
prob = osqp;

% Setup workspace
prob.setup(P, q, A, l, u, 'warm_start', true);

% solve it
tic,
out = prob.solve();
elapsedTime2 = toc;


% % Solve problem for different values of gamma parameter
% for i = 1 : length(gammas)
%     % Update linear cost
%     gamma = gammas(i);
%     q_new = [zeros(n+m,1); gamma*ones(n,1)];
%     prob.update('q', q_new);
% 
%     % Solve
%     res = prob.solve();
% end
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
