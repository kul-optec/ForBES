% demo_linprog
% demonstrates how to use ForBES to solve the linear program (LP)
%
% minimize c'*x subject to A*x=b, x>=0
%
% ForBES finds the minimum norm primal-dual pair by projecting the zero
% vector to the set defining the primal-dual optimality conditions for the
% LP

close all;
clear;

rng(0,'twister');

% Create linear program with known solution
m = 1000;
n = 1500;
A = sprandn(m, n, 100/n);
xopt = (rand(n, 1).*(rand(n, 1) < 0.5));
sopt = (rand(n, 1).*(xopt == 0));
yopt = rand(m, 1);
b = A*xopt;
c = A'*yopt + sopt;
costopt = c'*xopt;
% eliminate null rows
normrow = max(abs(A),[],2);
A = A(normrow>0,:);
b = b(normrow>0,1);
[m, n] = size(A);

%%                      run linprog
tic;[x,cost,exitflag,output,lambda]=linprog(c,[],[],A,b,zeros(n,1));timem = toc;
y = lambda.lower;

%%              run primal-dual LP ForBES solver
% dual problem is
%           maximize b'y subject to A'y+s = b, s >= 0
% Matrix K and vector d describe the affine subspace
% A*x = b, A'*y+s = 0, c'*x = b'*y
K = [sparse(n, n), A', speye(n); A, sparse(m, m+n); c', -b', sparse(1, n)];
d = [c; b; 0];
% setup forbes problem
% f1 is the sum of the squared norm plus the indicator of the affine
% subspace
prob.f1 = proj_over_affine(sparse(2*n+m,1),K,d);
% g is indicator of x>= 0, s>=0 
prob.g = indFreePos([1:n,n+m+1:2*n+m]);
prob.A1 = 1;
prob.B = -1;
prob.b = zeros(2*n+m,1);
% run forbes
tic;out = forbes(prob);timef = toc;
xf = out.x1(1:n);
yf = out.x1(n+1:n+m);
costf = c'*xf;
% run fast dual proximal gradient method (fast amm)
opt.fast = 1;opt.display = 0;
tic;outg = amm(prob,opt);timeg = toc;
xg = outg.x1(1:n);
yg = outg.x1(n+1:n+m);
costg = c'*xg;

fprintf('linprog  time elapsed (sec) : %.2f\n', timem);
fprintf('ForBES   time elapsed (sec) : %.2f\n', timef);
fprintf('Fast AMM time elapsed (sec) : %.2f\n\n', timeg);

fprintf('ForBES iterations   : %.0f\n', out.iterations);
fprintf('Fast AMM iterations : %.0f\n', outg.iterations);

