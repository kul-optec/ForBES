close all;
clear;

NUM_TOL_VAL = 1e-8;
NUM_TOL_DER = 1e-8;

A = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
b = [1; 2; 3];

[m, n] = size(A);

%% Lasso

f = forbes.functions.SqrNormL2();
lam = 1.0;
g = forbes.functions.NormL1(lam);
x0 = ones(n, 1);

prob = forbes.problems.ProblemComposite(f, A, -b, [], [], [], g, [], [], x0);
ops = forbes.fbe.FBOperations();

gams = [10.0/200, 5.0/200, 2.0/200, 1.0/200];

for ix = 1:10 % try several starting points

x = randn(n, 1);

for igam = 1:length(gams)

gam = gams(igam);
cache = forbes.fbe.FBCache(prob, x, gam, ops);

for idir = 1:10 % try several random directions

dir1 = randn(n, 1);
dir2 = -cache.Get_FPR();
cache.Set_Directions(dir1);
cache.Set_Directions([], dir2);

taus = [1.0, 0.5, 0.25, 0.125];

for itau = 1:length(taus)

    tau = taus(itau);

    cache_1 = cache.Get_CacheSegment(tau);
    cache_2 = forbes.fbe.FBCache(prob, x+tau*dir1+(1-tau)*dir2, gam, ops);

    assert(abs(cache_1.Get_FBE() - cache_2.Get_FBE())/abs(cache_2.Get_FBE()) <= NUM_TOL_VAL);

end

cache_1 = cache.Get_CacheSegment(0.0);
assert(norm(cache_1.Get_Point() - cache.Get_ProxGradStep(), inf) <= 1e-12);

end
end
end
