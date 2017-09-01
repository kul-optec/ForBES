close all;
clear;

ASSERT_EPS = 1e-10;
N_TESTS = 20;

% dense Q

n = 100;
Q = rand(n,n);
Q = Q'*Q + 1e-1*eye(n);
q = randn(n, 1);

f = forbes.functions.Quadratic(Q, q);
fc = forbes.functions.Conjugate(f);

for i = 1:N_TESTS
    x = randn(n, 1);
    y = randn(n, 1);
    [grad_f_x, f_x] = f.gradient(x);
    [grad_fc_y, fc_y] = fc.gradient(y);
    assert(f_x + fc_y >= x'*y);
    assert(abs(f_x + fc.call(grad_f_x) - x'*grad_f_x) <= ASSERT_EPS);
    assert(abs(fc_y + f.call(grad_fc_y) - grad_fc_y'*y) <= ASSERT_EPS);
end

% sparse Q

n = 500;
dens = 0.1;
rc = 0.1;
Q = sprandsym(n, dens, rc, 1) + 1e-1*speye(n);
q = randn(n, 1);

f = forbes.functions.Quadratic(Q, q);
fc = forbes.functions.Conjugate(f);

for i = 1:N_TESTS
    x = randn(n, 1);
    y = randn(n, 1);
    [grad_f_x, f_x] = f.gradient(x);
    [grad_fc_y, fc_y] = fc.gradient(y);
    assert(f_x + fc_y >= x'*y);
    assert(abs(f_x + fc.call(grad_f_x) - x'*grad_f_x) <= ASSERT_EPS);
    assert(abs(fc_y + f.call(grad_fc_y) - grad_fc_y'*y) <= ASSERT_EPS);
end
