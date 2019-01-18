close all;
clear;

myFolder = 'maros_meszaros_data';

if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

[Filename, pathname] = uigetfile('*.mat',myFolder,'MultiSelect','on');
matFiles             = fullfile(pathname,Filename);

ll       = length(matFiles);
out_maha = cell(1, ll);
out_osqp = cell(1, ll);
new      = {};

fprintf("##########################################\n");
for k = 1:ll
    baseFileName = Filename{k};
    new{k}       = char(baseFileName(1:end-4));
    fullFileName = fullfile(baseFileName);
    
    fprintf(1, '\n\nNow reading %s\n', fullFileName);
    matData{k} = load(fullFileName);
    dim        = matData{k}.n;
    P          = matData{k}.P + 1e-6*speye(dim);
    q          = matData{k}.q;
    lb         = matData{k}.l;
    ub         = matData{k}.u;
    A          = matData{k}.A;    
    y0         = zeros(matData{k}.m, 1);
   
    opt.memory  = 5;
    opt.maxit   = 4000;
    opt.tol     = 1e-3;
    opt.rel_tol = 1e-3;
    opt.abs_tol = 1e-3;
   
    opt.prescale   = 1;
    opt.adaptive   = 1;
    opt.solver     = 'nama';
    opt.linesearch = 'backtracking';
    opt.term       = @(prob, it, cache_0, cache_x, scale, Hessian) custom_termination(prob, it, cache_0, cache_x, scale,Hessian);
    out_maha{k}    = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt, y0);
    
    opt.solver     = 'zerofpr2';
    opt.linesearch = 'backtracking';
    opt.term       = @(prob, it, cache_0, cache_x, scale, Hessian) custom_termination(prob, it, cache_0, cache_x, scale,Hessian);
    out_panoc{k}   = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt, y0);
    
    
    fprintf("OSQP\n");
    prob      = osqp;
    fprintf("----- No Scaling, No Warm-start, No solution-polishing ------\n");
    prob.setup(P, q, A, lb, ub, 'rho',1.0,'alpha',1.6);
    out_osqp{k} = prob.solve();
     
    S(k,1).site_name = strcat('p-res','-', char(new{k}));
    S(k,2).site_name = strcat('Laug','-', char(new{k}));   
    S(k,1).value0  =  out_maha{k}.solver.solver.residual;
    S(k,2).value0  = -out_maha{k}.solver.solver.objective;
    S(k,1).value1  =  out_panoc{k}.solver.solver.residual;
    S(k,2).value1  = -out_panoc{k}.solver.solver.objective;
     
    row = 2;
    col = k;
    for kk = 1:numel(S)
        subplot(row,col,kk);
        semilogy(S(kk).value0);
        hold on
        semilogy(S(kk).value1);
        hold on
        legend('NAMA','P-NAMA')
        hold off
        title(S(kk).site_name)
    end
end

function s = printString(x)
    if floor(x) == x
        s = sprintf('%i', x);
    else
       s = sprintf('%0.2e', x); 
    end 
end  
   