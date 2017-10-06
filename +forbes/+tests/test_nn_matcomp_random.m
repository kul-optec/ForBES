rng(0);

m = 30; % number of rows
n = 30; % number of column of the original matrix M
d = 0.5; % density of coefficients sampled from M
r = 3; % rank of M

U = randn(m, r);
V = randn(n, r);
M = U*V';

P = sprand(m, n, d) ~= 0; % sampling pattern
B = full(M.*P);

f = forbes.functions.SqrNormL2(P);
aff = {1, -B};
lam = 2;
g = forbes.functions.NuclearNorm(lam, 'exact');
x0 = zeros(m, n);

TOL = 1e-8;

sol_FBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL));

x_star = sol_FBS.solution();

sol_AFBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL, 'fast', true));
assert(norm(sol_AFBS.solution - x_star, 'fro') <= 10*TOL*norm(x_star, 'fro'));

sol_NAMA = forbes(f, g, x0, aff, {}, forbes.solvers.NAMA('tol', TOL));
assert(norm(sol_NAMA.solution - x_star, 'fro') <= 10*TOL*norm(x_star, 'fro'));
