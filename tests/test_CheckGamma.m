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

Lf = norm(A)^2;

gam = 10/Lf;
bet = 0.05;

cache = FBCache(prob, x0, gam, ops);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 0);

gam = 0.9/Lf;

cache.Set_Gamma(gam);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 1);

%% Sparse logistic regression

f = logLoss(3.0);
lam = 10.0;
g = l1Norm(lam);
x0 = ones(n, 1);

prob = ProblemComposite([], [], [], f, A, -b, g, [], [], x0);

ops = FBOperations();

gam = 100.0/600;
bet = 0.05;

cache = FBCache(prob, x0, gam, ops);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 0);

gam = 1.0/600;

cache.Set_Gamma(gam);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 1);
