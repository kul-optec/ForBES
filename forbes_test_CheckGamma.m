function forbes_test_CheckGamma()

A = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
b = [1; 2; 3];

[m, n] = size(A);

f = quadLoss(1, zeros(m,1));
aff = {A, -b};
lam = 1.0;
g = l1Norm(lam);
x0 = zeros(n, 1);

opt.tol = 1e-6; % to create a structure

opt = ProcessOptions(opt);
prob = MakeProblem(f, g, x0, aff, [], opt);
prob = ProcessCompositeProblem(prob, opt);

x = ones(n, 1);
gam = 1.0/200;
bet = 0.05;

cache = CacheInit(prob, x, gam);
[flag, ~, ~] = CheckGamma(cache, gam, bet);

assert(flag == 1);
