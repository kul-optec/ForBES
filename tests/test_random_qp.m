%% makes a big problem
m = 500;
n = 1000;
A  = sparse(randn(m,n));
lA = -rand(m,1) * 2;
uA = +rand(m,1) * 2;
P = sprandsym(n,0.1,0.1,2) + speye(n);
q = randn(n,1);

%% QUADPROG SANITY CHECK

Aqp = [A;-A];
bqp = [uA;-lA];
tic
[x,f,flag] = quadprog(P,q,Aqp,bqp);
elapsedTime1 = toc;

%% OSQP SOLVER TEST
%make a solver object
osqpSolver = osqp;

%setup solver with data
osqpSolver.setup(P,q,A,lA,uA,'alpha', 1);

% solve it
tic,
out = osqpSolver.solve();
elapsedTime2 = toc;

%% NAMA Solver Test
opt.solver = 'nama';
opt.term = @(prob, it, cache_0, cache_x,scale) custom_termination(prob, it, cache_0, cache_x,scale);
tic,
nama_out = forbes_qp(P, q, A, lA, uA,[],[], [], [], opt);
elapsedTime3 = toc;


fprintf('\nThe time for the solvers\n');
fprintf('quadprog: %.2e, osqp:%.2e, nama:%.2e', elapsedTime1, elapsedTime2, elapsedTime3);
fprintf('osqp:(%.2e,%2d), nama:(%.2e,%2d)\n', out.info.pri_res, out.info.iter, nama_out.solver.solver.residual(end),nama_out.solver.solver.iterations);


return;




