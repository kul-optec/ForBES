A = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
b = [1; 2; 3];

[m, n] = size(A);

%% Lasso

f = quadLoss(1, zeros(m,1));
aff = {A, -b};
lam = 1.0;
g = l1Norm(lam);
x0 = ones(n, 1);

opt.maxit = 1;
opt.display = 0;
out = forbes(f, g, x0, aff, [], opt);
opt = out.opt;
prob = out.prob;

x = ones(n, 1);
gam = 10.0/200;
bet = 0.05;

cache = Cache_Init(prob, x, gam);
[flag, ~, ~] = Cache_CheckGamma(cache, gam, bet);

assert(flag == 0);

gam = 1.0/200;

cache = Cache_Init(prob, x, gam);
[flag, ~, ~] = Cache_CheckGamma(cache, gam, bet);

assert(flag == 1);

%% Sparse logistic regression

f = logLoss(3.0);
aff = {A, -b};
lam = 10.0;
g = l1Norm(lam);
x0 = ones(n, 1);

opt.maxit = 1;
opt.display = 0;
out = forbes(f, g, x0, aff, [], opt);
opt = out.opt;
prob = out.prob;

x = ones(n, 1);
gam = 100.0/600;
bet = 0.05;

cache = Cache_Init(prob, x, gam);
[flag, ~, ~] = Cache_CheckGamma(cache, gam, bet);

assert(flag == 0);

gam = 1.0/600;

cache = Cache_Init(prob, x, gam);
[flag, ~, ~] = Cache_CheckGamma(cache, gam, bet);

assert(flag == 1);
