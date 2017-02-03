function [x, y, s, info] = forbes_cp(data, K, opt)

info.name = 'FBCS matlab (ForBES wrapper)';
info.version = '0.1';

if nargin < 2
    x = []; y = []; s = [];
    return;
end

% some constants, flags etc.
normalize = 1;
scale = 1;
% verb_gap = 200;
% lbfgs_mem = 20;
% flag_small_stepsize = 0;
% eta = 0.00; % for nonmonotone line-search (eta = 0 makes it monotone, elsewhere we have eta = 0.85)
% max_iters_ls = 32;

% t0 = tic();

K = validate_cone(K);

n = length(data.c);
m = length(data.b);

unscaled_b = data.b;
unscaled_c = data.c;

unscaled_nm_b = norm(unscaled_b);
unscaled_nm_c = norm(unscaled_c);

data_orig = data;

work = struct();
if (normalize)
    [data, work] = normalize_data(data, K, scale, work);
    D = work.D;
    E = work.E;
    sc_b = work.sc_b;
    sc_c = work.sc_c;
else
    scale = 1;
    D = ones(m,1);
    E = ones(n,1);
    sc_b = 1;
    sc_c = 1;
end

% unwrap problem data
A = data.A;
b = data.b;
c = data.c;

% determine gamma
L = [sparse(n,n), A'; -A, sparse(m,m); -c', -b'];
eigsOpts.issym = 1;
eigsOpts.tol = 1e-3;
sqNormL = eigs(@(x)L*(L'*x), n+m+1, 1, 'LM', eigsOpts);

f = smoothCP(n, K, zeros(m+n,1), 1);
g = nonsmoothCP(m, n, K);
d = [-c; -b; 0];
opt.Lf = sqNormL;

t0 = tic();
out_forbes = forbes(f, g, zeros(n+m+1,1), [], {L, -1, d}, opt);
ttot = toc(t0);

% temporary
if opt.display > 0
    fprintf('Solver used : %s\n', out_forbes.name);
    fprintf('Iterations  : %d\n', out_forbes.iterations);
    fprintf('Stepsize    : %.2e (1/Lip = %.2e)\n', out_forbes.gam, 1/sqNormL);
    fprintf('Solve time  : %.2f\n', ttot);
end

x = out_forbes.x2(1:n);
y = out_forbes.x2(n+1:end);
s = out_forbes.z(n+1:n+m);

if (normalize)
    y = y ./ (D * sc_c);
    x = x ./ (E * sc_b);
    s = s .* (D / (sc_b * scale));
end

A_orig = data_orig.A;
b_orig = data_orig.b;
c_orig = data_orig.c;

unscaled_dres = A_orig'*y + c_orig;
unscaled_pres = b_orig - A_orig*x - s;
unscaled_cx = c_orig'*x;
unscaled_by = b_orig'*y;
unscaled_gap = -(unscaled_by + unscaled_cx);

info.resPri = norm(unscaled_pres)/(1 + unscaled_nm_b);
info.resDual = norm(unscaled_dres)/(1 + unscaled_nm_c);
info.relGap = abs(unscaled_gap)/(1 + abs(unscaled_cx) + abs(unscaled_by));

info.iter = out_forbes.iterations;
info.status = 'solved';

end
