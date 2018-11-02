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
for k = 1:ll
    baseFileName = Filename{k};
    new{k}       = string(baseFileName(1:end-4));
    fullFileName = fullfile(baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData{k}   = load(fullFileName);
    dim          = matData{k}.n;
    
   
    P          = matData{k}.P + 1e-3*speye(dim);
    q          = matData{k}.q;
    lb         = matData{k}.l;
    ub         = matData{k}.u;
    A          = matData{k}.A;
    opt.solver = 'nama';
   % opt.tol   = 1e-5;
    opt.term   = @(prob, it, cache_0, cache_x, scale) custom_termination(prob, it, cache_0, cache_x, scale); 
%   y0         = zeros(matData{k}.m, 1);
%   opt1       = opt;
%   opt1.maxit = 1;
%   inital_out = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt1);
    out_maha{k} = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt);
    
    
    %%%%  AMA
    opt.solver   = 'fbs';
    opt.method   = '';
    opt.adaptive = 1;
    opt.term     = @(prob, it, cache_0, cache_x, scale) custom_termination(prob, it, cache_0, cache_x, scale);
    out_maha1{k} = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt);
    
    
    err1     = norm(min(A*out_maha{k}.x-lb,0),'inf');
    err_nama = max(err1,norm(max(A*out_maha{k}.x-ub,0),'inf'));
    
    
    % Create an OSQP object
    prob1 = osqp;
    
    % Setup workspace and change alpha parameter
    prob1.setup(P, q, A, lb, ub, 'alpha', 1);
    
    % Solve problem
    out_osqp{k} = prob1.solve();%res is a struct and it have info structure
    
    err         = norm(min(A*out_osqp{k}.x-lb,0),'inf');
    err_osqp    = max(err,norm(max(A*out_osqp{k}.x-ub,0),'inf'));
   
    
%     %results- nonscaler struct array
% %     some = new{k};
%     results(k,1) = out_maha{k}.solver.solver.gam;
%     results(k,2) = out_maha{k}.solver.solver.time;
%     results(k,3) = out_osqp{k}.info.run_time;
%     results(k,4) = out_maha{k}.solver.solver.iterations;
%     results(k,5) = out_osqp{k}.info.iter;
%     results(k,6) = out_maha{k}.solver.solver.residual(end);
%     results(k,7) = out_osqp{k}.info.pri_res;
%     results(k,8) = -out_maha{k}.solver.solver.objective(end);
%     results(k,9) = out_osqp{k}.info.obj_val;
%     results(k,10) = out_maha{k}.status;
%     results(k,11) = out_osqp{k}.info.status_val;
%     results(k,12) = err_nama;
%     results(k,13) = err_osqp;
%     
%     csvwrite('results.csv',results)
  
    S(k,1).site_name = strcat('p-res','-', char(new{k}));
    S(k,2).site_name = strcat('obj','-', char(new{k}));   
    
    
    S(k,1).value = out_maha{k}.solver.solver.residual; 
    S(k,2).value = -out_maha{k}.solver.solver.objective;
    
    S(k,1).value1 =  out_maha1{k}.solver.solver.residual; 
    S(k,2).value1 = -out_maha1{k}.solver.solver.objective;
    
    
    row = 2;
    col = k;
    for kk = 1:numel(S)
        subplot(row,col,kk);
        loglog(S(kk).value);
        hold on 
        loglog(S(kk).value1);
        hold off
        title(S(kk).site_name)
    end
    
end
