close all;
clear;

rng(0, 'twister');

lam_fact = 0.05;

%% load from file

PATH_TO_DATASETS = '/Users/lorenzo/Datasets/L1_Testset_mat/';
problem = 1;
shortname = strcat('spear_inst_', num2str(problem));
filename = strcat(PATH_TO_DATASETS, shortname, '.mat');
data = load(filename, '-mat');
A = data.A;
b = data.b;
[m, n] = size(A);
lam_max = norm(A'*b,'inf');
lam = lam_fact*lam_max;

%% set up problem

f = quadLoss(1, zeros(m, 1));
g = l1Norm(lam);
x0 = zeros(n, 1);

%% set stopping criterion and recording options

eigsOpt.issym = 1;
eigsOpt.tol = 1e-3;
funHessian = @(x) A'*(A*x);
Lf = eigs(funHessian, n, 1, 'LM', eigsOpt);

opt.Lf = Lf;
opt.maxit = 2000;
opt.memory = 5;
opt.display = 1;
opt.tol = 1e-8;

%% define algorithms to run

opts = {};
opts{end+1} = opt; opts{end}.solver = 'fbs'; opts{end}.variant = 'fast';
opts{end+1} = opt; opts{end}.solver = 'minfbe'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = opt; opts{end}.solver = 'zerofpr'; opts{end}.method = 'lbfgs'; opts{end}.linesearch = 'backtracking';
opts{end+1} = opt; opts{end}.solver = 'zerofpr_new'; opts{end}.method = 'lbfgs';

%% run algorithms

outs = {};
names = {};
for k=1:length(opts)
    outs{k} = forbes(f, g, x0, {A, -b}, [], opts{k});
    names{k} = outs{k}.name;
    fprintf('%30s, iters: %5d, matvecs: %5d, prox: %5d, time: %.3f\n', ...
      names{k}, outs{k}.iterations, outs{k}.operations.C1, outs{k}.operations.proxg, outs{k}.ts(end));
end
