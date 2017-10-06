rng(0);

A = [1,  2, -1, -1; ...
    -2, -1,  0, -1; ...
    3,  0,  4, -1; ...
    -4, -1, -3,  1; ...
    5,  3,  2,  3]';
b = [1, 2, 3, 4]';

[m, n] = size(A);

f = forbes.functions.SqrNormL2();
aff = {A, -b};
lam = 0.1*norm(A'*b, 'inf');
g = forbes.functions.NormL1(lam);
x0 = zeros(n, 1);

x_star = [-3.877278911564627e-01; 0; 0; 2.174149659863943e-02; 6.168435374149660e-01];

TOL = 1e-8;

sol_FBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL));
assert(norm(sol_FBS.solution() - x_star, 'inf') <= 10*TOL*norm(x_star, 'inf'));

sol_AFBS = forbes(f, g, x0, aff, {}, forbes.solvers.FBS('tol', TOL, 'fast', true));
assert(norm(sol_AFBS.solution() - x_star, 'inf') <= 10*TOL*norm(x_star, 'inf'));

sol_NAMA = forbes(f, g, x0, aff, {}, forbes.solvers.NAMA('tol', TOL));
assert(norm(sol_NAMA.solution() - x_star, 'inf') <= 10*TOL*norm(x_star, 'inf'));
