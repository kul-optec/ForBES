rng(0);

B = % TODO: put something here
P = B ~= 0;

f = forbes.functions.SqrNormL2(P);
aff = {1, -B};
lam = % TODO: put something here
g = forbes.functions.NuclearNorm(lam, 'exact');
x0 = zeros(m, n);

x_star = % TODO: put something here

TOL = 1e-8;

sol_FBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL));
assert(norm(sol_FBS.solution() - x_star, 'fro') <= 10*TOL*norm(x_star, 'fro'));

sol_AFBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL, 'fast', true));
assert(norm(sol_AFBS.solution() - x_star, 'fro') <= 10*TOL*norm(x_star, 'fro'));

sol_NAMA = forbes(f, g, x0, aff, {}, forbes.solvers.NAMA('tol', TOL));
assert(norm(sol_NAMA.solution() - x_star, 'fro') <= 10*TOL*norm(x_star, 'fro'));
