NUM_TOL_VAL = 1e-8;
NUM_TOL_DER = 1e-8;

A = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
b = [1; 2; 3];

[m, n] = size(A);

%% Lasso

f = quadLoss(1, zeros(m,1));
aff = {A, -b};
lam = 1.0;
g = l1Norm(lam);
x0 = ones(n, 1);

opt.tol = 1e-6; % to create a structure

opt = Process_Options(opt);
prob = Process_MakeProblem(f, g, x0, aff, []);
prob = Process_CompositeProblem(prob, opt);

gams = [10.0/200, 5.0/200, 2.0/200, 1.0/200];

for ix = 1:10 % try several starting points

x = randn(n, 1);

for igam = 1:length(gams)

gam = gams(igam);
cache = Cache_Init(prob, x, gam);

for idir = 1:10 % try several random directions

dir = randn(n, 1);
cache = Cache_LineSearch(cache, dir);

taus = [1.0, 0.5, 0.25, 0.125];

for itau = 1:length(taus)

    tau = taus(itau);

    cache_1 = Cache_LineFBE(cache, tau, 1);
    cache_2 = Cache_Init(prob, x+tau*dir, gam);
    cache_2 = Cache_FBE(cache_2, gam);

    assert(abs(cache_1.FBE - cache_2.FBE)/abs(cache_2.FBE) <= NUM_TOL_VAL);

    cache_1 = Cache_LineFBE(cache, tau, 2);
    cache_2 = Cache_GradFBE(cache_2, gam);
    slope = cache_2.gradFBE'*dir;

    assert(abs(cache_1.dFBE - slope)/abs(slope) <= NUM_TOL_DER);

    cache_1 = Cache_LineFBE(cache, tau, 3);

    assert(abs(cache_1.FBE - cache_2.FBE)/abs(cache_2.FBE) <= NUM_TOL_VAL);
    assert(abs(cache_1.dFBE - slope)/abs(slope) <= NUM_TOL_DER);

end
end
end
end
