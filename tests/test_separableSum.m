close all;
clear;

ASSERT_EPS = 1e-14;

% two simple functions

w = rand(10,1);
f1 = forbes.functions.SqrNormL2(w); % R^10 -> R

mu = 1;
f2 = forbes.functions.LogisticLoss(mu); % R^? -> R

fSum = forbes.functions.SeparableSum({f1, f2}, {10, 20}); % R^30 -> R

for i = 1:100
    x = randn(30, 1);
    [grad1, v1] = fSum.gradient(x);
    v2 = 0.5*((w.*(x(1:10)))'*(x(1:10))) + mu*sum(log(1+exp(-x(11:end))));
    grad2 = [vec(w.*(x(1:10))); -mu*exp(-x(11:end))./(1+exp(-x(11:end)))];
    assert(abs(v1-v2)/(1+abs(v2)) <= ASSERT_EPS);
    assert(norm(grad1-grad2, inf)/(1+norm(grad2, inf)) <= ASSERT_EPS);
end

% two less simple functions

f1 = forbes.functions.SqrNormL2(1); % R^? -> R

mu = 1;
f2 = forbes.functions.LogisticLoss(mu); % R^? -> R

fSum = forbes.functions.SeparableSum({f1, f2}, {[20, 30], 50}); % R^650 -> R

for i = 1:100
    x = randn(650, 1);
    [grad1, v1] = fSum.gradient(x);
    v2 = 0.5*norm(reshape(x(1:600), 20, 30),'fro')^2 + mu*sum(log(1+exp(-x(601:end))));
    grad2 = [vec(reshape(x(1:600), 20, 30)); -mu*exp(-x(601:end))./(1+exp(-x(601:end)))];
    assert(norm(v1-v2, inf)/(1+abs(v2)) <= ASSERT_EPS);
    assert(norm(grad1-grad2, inf)/(1+norm(grad2, inf)) <= ASSERT_EPS);
end

% composition with affine mappings
