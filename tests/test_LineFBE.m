close all;
clear;

NUM_TOL_VAL = 1e-8;
NUM_TOL_DER = 1e-8;

A = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
b = [1; 2; 3];

[m, n] = size(A);

%% Lasso

f = quadLoss(1, zeros(m,1));
lam = 1.0;
g = l1Norm(lam);
x0 = ones(n, 1);

prob = ProblemComposite(f, A, -b, [], [], [], g, [], [], x0);
ops = FBOperations();

gams = [10.0/200, 5.0/200, 2.0/200, 1.0/200];

for ix = 1:10 % try several starting points

x = randn(n, 1);

for igam = 1:length(gams)

gam = gams(igam);
cache = FBCache(prob, x, gam, ops);

for idir = 1:10 % try several random directions

dir = randn(n, 1);
cache.Set_Directions(dir);

taus = [1.0, 0.5, 0.25, 0.125];

for itau = 1:length(taus)

    tau = taus(itau);

    cache_1 = cache.Get_CacheLine(tau, 1);
    cache_2 = FBCache(prob, x+tau*dir, gam, ops);

    assert(abs(cache_1.Get_FBE() - cache_2.Get_FBE())/abs(cache_2.Get_FBE()) <= NUM_TOL_VAL);

    cache_1 = cache.Get_CacheLine(tau, 2);
    gradFBE = cache_2.Get_GradFBE();
    slope = gradFBE'*dir;

    assert(abs(cache_1.dFBE - slope)/abs(slope) <= NUM_TOL_DER);

    cache_1 = cache.Get_CacheLine(tau, 3);

    assert(abs(cache_1.Get_FBE() - cache_2.Get_FBE())/abs(cache_2.Get_FBE()) <= NUM_TOL_VAL);
    assert(abs(cache_1.Get_Slope() - slope)/abs(slope) <= NUM_TOL_DER);

end
end
end
end
