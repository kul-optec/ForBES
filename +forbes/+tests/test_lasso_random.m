rng(0);

m = 50;
n = 200;

A = randn(m, n);
b = randn(m, 1);

f = forbes.functions.SqrNormL2();
aff = {A, -b};
lam = 0.3*norm(A'*b, 'inf');
g = forbes.functions.NormL1(lam);
x0 = zeros(n, 1);

TOL = 1e-8;

sol_FBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL));

x_star = sol_FBS.solution();

sol_AFBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL, 'fast', true));
assert(norm(sol_AFBS.solution() - x_star, 'inf') <= 10*TOL*norm(x_star, 'inf'));

sol_NAMA = forbes(f, g, x0, aff, {}, forbes.solvers.NAMA('tol', TOL));
assert(norm(sol_NAMA.solution() - x_star, 'inf') <= 10*TOL*norm(x_star, 'inf'));
