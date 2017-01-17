close all;
clear;

A = randn(15,10);
Q = A*A';
q = randn(15, 1);
f = quadratic(Q, q);
g = indBox(-1, 1);
x0 = randn(15,1);
opt.tol = 1e-8;
opt.display = 0;
out = forbes(f, g, x0, [], [], opt);

fs = {quadLoss(), logLoss(), quadLoss()};
gs = {indPos(), l1Norm()};
aff = { randn(5, 3), randn(5, 6), zeros(5, 1);
        randn(7, 3), randn(7, 6), randn(7, 1);
        randn(9, 3), randn(9, 6), -ones(9, 1);};
x0 = randn(9,1);
opt.tol = 1e-8;
opt.display = 0;
out = forbes(fs, gs, x0, aff, [], opt);
