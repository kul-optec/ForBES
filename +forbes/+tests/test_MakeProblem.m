close all;
clear;

A = randn(15,10);
Q = A*A';
q = randn(15, 1);
f = forbes.functions.Quadratic(Q, q);
g = forbes.functions.IndBox(-1, 1);
x0 = randn(15,1);
opt.tol = 1e-8;
opt.display = 0;
out = forbes(f, g, x0, [], [], opt);

fs = {forbes.functions.SqrNormL2(), forbes.functions.LogisticLoss(), forbes.functions.SqrNormL2()};
gs = {forbes.functions.IndPos(), forbes.functions.NormL1()};
aff = { randn(5, 3), randn(5, 6), zeros(5, 1);
        randn(7, 3), randn(7, 6), randn(7, 1);
        randn(9, 3), randn(9, 6), -ones(9, 1);};
x0 = randn(9,1);
opt.tol = 1e-8;
opt.display = 0;
out = forbes(fs, gs, x0, aff, [], opt);
