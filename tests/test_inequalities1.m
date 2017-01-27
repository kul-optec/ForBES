close all;
clear;

N_TESTS = 2000;

m = 50;
n = 200;
A = randn(m, n);
b = randn(m, 1);
f = quadLoss();
aff = {A, -b};
lam = 0.3*norm(A'*b, 'inf');
g = distBox(-1, 1, lam);
x0 = zeros(n, 1);

prob = ProblemComposite(f, A, -b, [], [], [], g, [], [], x0);
opt.adaptive = 0;
opt.beta = 0.05;
[Lf, adaptive] = prob.Get_Lipschitz(opt);
gam = (1-opt.beta)/Lf;
ops = FBOperations();

for i = 1:N_TESTS
  x = 3*randn(n, 1);
  cache_x = FBCache(prob, x, gam, ops);
  z = cache_x.Get_ProxGradStep();
  cache_z = FBCache(prob, z, gam, ops);
  assert(cache_z.Get_f() + cache_x.gz <= cache_x.Get_FBE() - opt.beta/(2*gam)*cache_x.Get_NormFPR()^2 + abs(cache_x.Get_FBE())*1e-12);
  assert(cache_z.Get_FBE() <= cache_z.Get_f() + cache_x.gz - 1/(2*gam)*cache_z.Get_NormFPR()^2 + abs(cache_z.Get_FBE())*1e-12);
end
