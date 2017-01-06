ASSERT_EPS = 1e-10;
N_TESTS = 20;

% dense Q

n = 100;
Q = rand(n,n);
Q = Q'*Q + 1e-1*eye(n);
q = randn(n, 1);

f = quadratic(Q, q);

call_f = f.makef();
call_fc = f.makefconj();

for i = 1:N_TESTS
    x = randn(n, 1);
    y = randn(n, 1);
    [f_x, grad_f_x] = call_f(x);
    [fc_y, grad_fc_y] = call_fc(y);
    assert(f_x + fc_y >= x'*y);
    assert(abs(f_x + call_fc(grad_f_x) - x'*grad_f_x) <= ASSERT_EPS);
    assert(abs(fc_y + call_f(grad_fc_y) - grad_fc_y'*y) <= ASSERT_EPS);
end

% sparse Q

n = 500;
dens = 0.1;
rc = 0.1;
Q = sprandsym(n, dens, rc, 1) + 1e-1*speye(n);
q = randn(n, 1);

f = quadratic(Q, q);

call_f = f.makef();
call_fc = f.makefconj();

for i = 1:N_TESTS
    x = randn(n, 1);
    y = randn(n, 1);
    [f_x, grad_f_x] = call_f(x);
    [fc_y, grad_fc_y] = call_fc(y);
    assert(f_x + fc_y >= x'*y);
    assert(abs(f_x + call_fc(grad_f_x) - x'*grad_f_x) <= ASSERT_EPS);
    assert(abs(fc_y + call_f(grad_fc_y) - grad_fc_y'*y) <= ASSERT_EPS);
end