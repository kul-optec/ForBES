close all;
clear;

N_TESTS = 2000;

% dense Q

n = 10;
m = 20;
Q = randn(n,n);
Q = 0.5*(Q'+Q);
Q = Q + (1e-1 - min(eig(Q)))*eye(n);
q = randn(n, 1);
A = randn(m, n);
low = -0.2; hi = 0.2;
f = forbes.functions.Quadratic(Q, q);
g = forbes.functions.NormL1();

y0 = randn(m, 1);

prob = forbes.problems.ProblemComposite(forbes.functions.Conjugate(f), -A', [], [], [], [], forbes.functions.Conjugate(g), 1, [], y0);
opt.adaptive = 0;
opt.beta = 0.05;
[Lf, adaptive] = prob.Get_Lipschitz(opt);
gam = (1-opt.beta)/Lf;
ops = forbes.fbe.FBOperations();

for i = 1:N_TESTS
  x = randn(m, 1);
  cache_x = forbes.fbe.FBCache(prob, x, gam, ops);
  z = cache_x.Get_ProxGradStep();
  cache_z = forbes.fbe.FBCache(prob, z, gam, ops);
  assert(cache_z.Get_f() + cache_x.gz <= cache_x.Get_FBE() - opt.beta/(2*gam)*cache_x.Get_NormFPR()^2 + abs(cache_x.Get_FBE())*1e-12);
  assert(cache_z.Get_FBE() <= cache_z.Get_f() + cache_x.gz - 1/(2*gam)*cache_z.Get_NormFPR()^2 + abs(cache_z.Get_FBE())*1e-12);
end
