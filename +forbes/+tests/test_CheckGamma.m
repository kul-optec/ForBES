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

Lf = norm(A)^2;

gam = 10/Lf;
bet = 0.05;

cache = forbes.fbe.FBCache(prob, x0, gam, ops);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 0);

gam = 0.9/Lf;

cache.Set_Gamma(gam);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 1);

%% Sparse logistic regression

f = forbes.functions.LogisticLoss(3.0);
lam = 10.0;
g = forbes.functions.NormL1(lam);
x0 = ones(n, 1);

prob = forbes.problems.ProblemComposite([], [], [], f, A, -b, g, [], [], x0);

ops = forbes.fbe.FBOperations();

gam = 100.0/600;
bet = 0.05;

cache = forbes.fbe.FBCache(prob, x0, gam, ops);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 0);

gam = 1.0/600;

cache.Set_Gamma(gam);
[flag, ~] = cache.Check_Gamma(bet);

assert(flag == 1);
