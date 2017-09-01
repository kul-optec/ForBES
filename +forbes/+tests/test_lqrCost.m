close all;
clear;

ASSERT_EPS = 1e-10;
N_TESTS = 20;

% Generate problem

n_x = 6;
n_u = 2;

A = randn(n_x, n_x); A = 2*(A/norm(A));
B = randn(n_x, n_u); B = B/norm(B);
Q = 1*eye(n_x);
R = 1e-1*eye(n_u);
Q_N = 1*Q;
N = 10;
x0 = randn(n_x, 1);

% Build up big constraint/cost matrices, so as to be able to
% compute f conjugate and its gradient via some QP solver
%
%   f(x) = (1/2)(x'*H*x) + ind{A*x = b}
%
%   fc(y) = sup_x { y'*x - f(x) }
%   => minimize_x f(x) - y'*x s.t. A*x = b
%   => obtain x_sol = grad_fc_y
%   => compute fc(y) = y'*x_sol - f(x_sol)

block_eq = [A, B, -eye(n_x)];
A_eq = sparse((N+1)*n_x, (N+1)*n_x + N*n_u);
diag_H = {};
for i = 0:N-1 % build up constraints
    basei = i*n_x;
    basej = i*(n_x+n_u);
    A_eq(basei+1:basei+n_x, basej+1:basej+2*n_x+n_u) = block_eq;
    diag_H{2*i+1} = Q;
    diag_H{2*i+2} = R;
end
diag_H{2*N+1} = Q_N;
H = blkdiag(diag_H{:});
b_eq = zeros((N+1)*n_x, 1);

% Test with no reference (zero reference)

f = forbes.functions.LQRCost(x0, Q, R, Q_N, A, B, N);
fc = forbes.functions.Conjugate(f);

for i=1:N_TESTS
    y = randn(N*(n_x+n_u)+n_x, 1);
    [grad_fc_y, fc_y] = fc.gradient(y);
    % test conjugate subgradient theorem
    fx = 0.5*(grad_fc_y'*H*grad_fc_y);
    assert(abs(grad_fc_y'*y - fx - fc_y) <= 1e-12*(1+abs(fc_y)));
    % evaluate gradient numerically
    grad_fc_y_num = forbes.utils.numdiff(@(x) fc.call(x), y);
    assert(norm(grad_fc_y-grad_fc_y_num, 'inf') <= 1e-6*(1+norm(grad_fc_y, 'inf')));
    % evaluate by solving a QP
%     opt_qp = optimoptions('quadprog','Display','off');
%     [grad_fc_y_qp] = quadprog(H,-y,[],[],A_eq,b_eq,[],[],[],opt_qp);
%     fc_y_qp = y'*grad_fc_y_qp - 0.5*(grad_fc_y_qp'*H*grad_fc_y_qp);
    % test equivalence
%     assert(abs(fc_y-fc_y_qp) <= 1e-8);
%     assert(norm(grad_fc_y-grad_fc_y_qp, 'inf') <= 1e-8);
end

% Test with reference state

xref = randn(n_x, 1);
tran = [repmat([Q*xref; zeros(n_u, 1)], N, 1); Q_N*xref];
f = forbes.functions.LQRCost(x0, Q, R, Q_N, A, B, N, xref);
fc = forbes.functions.Conjugate(f);

for i=1:N_TESTS
    y = randn(N*(n_x+n_u)+n_x, 1);
    [grad_fc_y, fc_y] = fc.gradient(y);
    % test conjugate subgradient theorem
    fx = 0.5*((grad_fc_y-tran)'*H*(grad_fc_y-tran));
    assert(abs(grad_fc_y'*y - fx - fc_y) <= 1e-12*(1+abs(fc_y)));
    % evaluate gradient numerically
    grad_fc_y_num = forbes.utils.numdiff(@(x) fc.call(x), y);
    assert(norm(grad_fc_y-grad_fc_y_num, 'inf') <= 1e-6*(1+norm(grad_fc_y, 'inf')));
    % evaluate by solving a QP
%     opt_qp = optimoptions('quadprog','Display','off');
%     [grad_fc_y_qp] = quadprog(H,-y,[],[],A_eq,b_eq,[],[],[],opt_qp);
%     fc_y_qp = y'*grad_fc_y_qp - 0.5*(grad_fc_y_qp'*H*grad_fc_y_qp);
    % test equivalence
%     assert(abs(fc_y-fc_y_qp) <= 1e-8);
%     assert(norm(grad_fc_y-grad_fc_y_qp, 'inf') <= 1e-8);
end
