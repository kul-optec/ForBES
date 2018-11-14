close all;
clear;

 myFolder = 'maros_meszaros_data';

if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end



[Filename, pathname] = uigetfile('*.mat',myFolder,'MultiSelect','on');
matFiles = fullfile(pathname,Filename);


ll = length(matFiles);
out_maha = cell(1, ll);
out_osqp = cell(1, ll);
new = {};
for k = 1:ll
    baseFileName = Filename{k};
    new{k} = string(baseFileName(1:end-4));
    
    fullFileName = fullfile(baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData{k} = load(fullFileName);
    dim = matData{k}.n;
    
    P  = matData{k}.P + 1e-3*speye(dim);
    q  = matData{k}.q;
    lb = matData{k}.l;
    ub = matData{k}.u;
    A  = matData{k}.A;
    opt.solver = 'nama';
    opt.memory = 10;
    opt.maxit = 10000;
    fprintf("-----No Scaling No adaptive ForBes-NAMA------\n");
    opt.prescale = false;
    opt.adaptive = false;
    out_maha2 = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt);
    errp_forbes_qp0 = max(  norm(min(A*out_maha2.x-lb,0),'inf'),norm(max(A*out_maha2.x-ub,0),'inf'));
    errd_forbes_qp0 = norm(P*out_maha2.x+q+A'*out_maha2.y_ineq,inf);
    errs_forbes_qp0 = max(norm(min(max(out_maha2.y_ineq,0),abs(ub-A*out_maha2.x)),inf),norm(min(-min(out_maha2.y_ineq,0),abs(A*out_maha2.x-lb)),inf));
    iter_forbes_qp0 = out_maha2.solver.solver.iterations;
    
    
    fprintf("-----Jacobi Scaling and adaptive ForBes-NAMA------\n");
    opt.prescale = true;% Jacobi
    opt.adaptive  = true;
    opt.term = @(prob, it, cache_0, cache_x, scale) custom_termination(prob, it, cache_0, cache_x, scale);
    out_maha{k} = forbes_qp(P, q, A, lb, ub,[],[], [], [], opt);
    errp_forbes_qp1 = max( norm(min(A*out_maha{k}.x-lb,0),'inf'),norm(max(A*out_maha{k}.x-ub,0),'inf'));
    errd_forbes_qp1 = norm(P*out_maha{k}.x+q+A'*out_maha{k}.y_ineq,inf);
    errs_forbes_qp1 = max(norm(min(max(out_maha{k}.y_ineq,0),abs(ub-A*out_maha{k}.x)),inf),norm(min(-min(out_maha{k}.y_ineq,0),abs(A*out_maha{k}.x-lb)),inf));
    iter_forbes_qp1 = out_maha{k}.solver.solver.iterations;
 
    
    fprintf("NAMA-QP\n");
    
    yinit = zeros(matData{k}.m, 1);
    opt.maxit = 10000;
    opt.memory = 10;
    opt.prescale = 0;% no scaling
    opt.adaptive = 0;
    [x00,z00,y00,iter00,res_pri00] = nama_qp(P,q,A,lb,ub,opt,yinit);
    errp_nama_qp00 = max(norm(min(A*x00-lb,0),'inf'),norm(max(A*x00-ub,0),'inf'));
    errd_nama_qp00 = norm(P*x00+q+A'*y00,inf);
    errs_nama_qp00 = max(norm(min(max(y00,0),abs(ub-A*x00)),inf),norm(min(-min(y00,0),abs(A*x00-lb)),inf));
    fprintf("No Scaling-no adaptive NAMA-QP %d, %.2e\n",iter00, res_pri00(end));
    
    
    
    yinit = zeros(matData{k}.m, 1);
    opt.maxit = 10000;
    opt.prescale = 0;% no scaling
    opt.adaptive = 1;
    [x0,z0,y0,iter0,res_pri0] = nama_qp(P,q,A,lb,ub,opt,yinit);
    errp_nama_qp0 = max(norm(min(A*x0-lb,0),'inf'),norm(max(A*x0-ub,0),'inf'));
    errd_nama_qp0 = norm(P*x0+q+A'*y0,inf);
    errs_nama_qp0 = max(norm(min(max(y0,0),abs(ub-A*x0)),inf),norm(min(-min(y0,0),abs(A*x0-lb)),inf));
    fprintf("No Scaling-but adaptive NAMA-QP %d, %.2e\n",iter0,res_pri0(end));
    
    
%     opt.prescale = 1;
%     opt.adaptive = 0;
%     [~,~,~,iter10,~,res_pri10] = nama_qp(P,q,A,lb,ub,opt,yinit);
%     fprintf("Jacobi Scaling-no adaptive NAMA-QP %d, %.2e\n",iter10,res_pri10(end));
%     
%     
%     opt.prescale = 1;% Jacobi
%     opt.adaptive = 1;
%     [x1,z1,y1,iter1,Laug1,res_pri1] = nama_qp(P,q,A,lb,ub,opt,yinit);
%     err_nama_qp1 = max( norm(min(A*x1-lb,0),'inf'),norm(max(A*x1-ub,0),'inf'));
%     errd_nama_qp1 = norm(P*x1+q+A'*y1,inf);
%     errs_nama_qp1 = max(norm(min(max(y1,0),abs(ub-A*x1)),inf),norm(min(-min(y1,0),abs(A*x1-lb)),inf));
%     fprintf("Jacobi Scaling and adaptive NAMA-QP %d, %.2e\n",iter1,res_pri1(end));
%     
%     
%     opt.prescale = 2;% Ruiz-non-adaptive
%     opt.adaptive = 0;
%     [x20,z20,y20,iter20,res_pri20] = nama_qp(P,q,A,lb,ub,opt,yinit);
%     err_nama_qp20  = max( norm(min(A*x20-lb,0),'inf'),norm(max(A*x20-ub,0),'inf'));
%     errd_nama_qp20 = norm(P*x20+q+A'*y20,inf);
%     errs_nama_qp20 = max(norm(min(max(y20,0),abs(ub-A*x20)),inf),norm(min(-min(y20,0),abs(A*x20-lb)),inf));
%     fprintf("Ruiz Scaling and no adaptive NAMA-QP %d, %.2e\n",iter20,res_pri20(end));
%     
%     
%     opt.prescale = 2;% Ruiz
%     opt.adaptive = 1;
%     [x2,z2,y2,iter2,res_pri2] = nama_qp(P,q,A,lb,ub,opt,yinit);
%     err_nama_qp2  = max( norm(min(A*x2-lb,0),'inf'),norm(max(A*x2-ub,0),'inf'));
%     errd_nama_qp2 = norm(P*x2+q+A'*y2,inf);
%     errs_nama_qp2 = max(norm(min(max(y2,0),abs(ub-A*x2)),inf),norm(min(-min(y2,0),abs(A*x2-lb)),inf));
%     fprintf("Ruiz Scaling and adaptive NAMA-QP %d, %.2e\n",iter2,res_pri2(end));
%    
%    
%     opt.prescale = 3;% Dual Ruiz-non-adaptive
%     opt.adaptive = 0;
%     [x30,z30,y30,iter30,res_pri30] = nama_qp(P,q,A,lb,ub,opt,yinit);
%     err_nama_qp30 = max( norm(min(A*x30-lb,0),'inf'),norm(max(A*x30-ub,0),'inf'));
%     errd_nama_qp30 = norm(P*x30+q+A'*y30,inf);
%     errs_nama_qp30 = max(norm(min(max(y30,0),abs(ub-A*x30)),inf),norm(min(-min(y30,0),abs(A*x30-lb)),inf));
%     fprintf("Dual Ruiz Scaling and no adaptive NAMA-QP %d, %.2e\n",iter30,res_pri30(end));
%     
%     
%     opt.prescale = 3;% Dual Ruiz
%     opt.adaptive = 1;
%     [x3,z3,y3,iter3,res_pri3] = nama_qp(P,q,A,lb,ub,opt,yinit);
%     err_nama_qp3 = max( norm(min(A*x3-lb,0),'inf'),norm(max(A*x3-ub,0),'inf'));
%     errd_nama_qp3 = norm(P*x3+q+A'*y3,inf);
%     errs_nama_qp3 = max(norm(min(max(y3,0),abs(ub-A*x3)),inf),norm(min(-min(y3,0),abs(A*x3-lb)),inf));
%     fprintf("Dual Ruiz Scaling and adaptive NAMA-QP %d, %.2e\n",iter3,res_pri3(end));
%     
% %     % Create an OSQP object
% %     prob0 = osqp;
% %     % Setup workspace and change alpha parameter
% %     fprintf("----- NO Scaling OSQP------\n");
% %     prob0.setup(P, q, A, lb, ub, 'rho',1.0,'alpha',1.6);
% %     
% %     % Solve problem
% %     out_osqp0{k} = prob0.solve();%res is a struct and it have info structure
% %     err = norm(min(A*out_osqp0{k}.x-lb,0),'inf');
% %     errp_osqp0 = max(err,norm(max(A*out_osqp0{k}.x-ub,0),'inf'));
% %     errd_osqp0 = norm(P*out_osqp0{k}.x+q+A'*out_osqp0{k}.y,inf);
% %     errs_osqp0 = max(norm(min(max(out_osqp0{k}.y,0),abs(ub-A*out_osqp0{k}.x)),inf),norm(min(-min(out_osqp0{k}.y,0),abs(A*out_osqp0{k}.x-lb)),inf));
% % 
%     % Create an OSQP object
%     prob1 = osqp;
%     % Setup workspace and change alpha parameter
%     fprintf("----- Ruiz Scaling OSQP------\n");
%     prob1.setup(P, q, A, lb, ub, 'rho',1.0,'alpha',1.6);
%     
%     % Solve problem
%     out_osqp{k} = prob1.solve();%res is a struct and it have info structure
%     err = norm(min(A*out_osqp{k}.x-lb,0),'inf');
%     errp_osqp = max(err,norm(max(A*out_osqp{k}.x-ub,0),'inf'));
%     errd_osqp = norm(P*out_osqp{k}.x+q+A'*out_osqp{k}.y,inf);
%     errs_osqp = max(norm(min(max(out_osqp{k}.y,0),abs(ub-A*out_osqp{k}.x)),inf),norm(min(-min(out_osqp{k}.y,0),abs(A*out_osqp{k}.x-lb)),inf));
%     
%    
%     results(k,1)  = out_maha{k}.solver.solver.gam;
%     results(k,2)  = out_maha{k}.solver.solver.iterations;
%     results(k,3)  = out_osqp{k}.info.iter;
%     results(k,4)  = iter00;
%     results(k,5)  = iter10;
%     results(k,6)  = iter1;
%     results(k,7)  = out_maha{k}.solver.solver.residual(end);
%     results(k,8)  = out_osqp{k}.info.pri_res;
%     results(k,9)  = res_pri00(end);
%     results(k,10) = res_pri1(end);
%     results(k,11) = res_pri11(end);
%     
%     
%     
%     %results- nonscaler struct array
% %     results(k,1)  = out_maha{k}.solver.solver.gam;
% %     results(k,2)  = out_maha{k}.solver.solver.time;
% %     results(k,3)  = out_osqp{k}.info.run_time;
% %     results(k,4)  = out_maha{k}.solver.solver.iterations;
% %     results(k,5)  = out_osqp{k}.info.iter;
% %     results(k,6)  = out_maha{k}.solver.solver.residual(end);
% %     results(k,7)  = out_osqp{k}.info.pri_res;
% %     results(k,8)  = -out_maha{k}.solver.solver.objective(end);
% %     results(k,9)  = out_osqp{k}.info.obj_val;
% %     results(k,10) = out_maha{k}.status;
% %     results(k,11) = out_osqp{k}.info.status_val;
% %     results(k,12) = errp_forbes_qp1;
% %     results(k,13) = errp_osqp;
%     Names = {'gam_j','iter_j','iter_o','iter_n','iter_nk','iter_nj','res_j','res_oo','res_jo','res_oq','res_o'};
%    % Names = {'gam_j','gam_r','iter_j','iter_r','iter_o','res_j', 'res_r','res_o','obj_j','obj_r','obj_o'};
%     T = cell2table(num2cell(results),'VariableNames',Names);%'RowNames',new,
%     writetable(T,'results.csv','WriteRowNames',true);
% 
% 
% 
%     S(k,1).site_name = strcat('p-res','-', char(new{k}));
%    % S(k,2).site_name = strcat('obj','-', char(new{k}));   
%     S(k,1).value = out_maha{k}.solver.solver.residual; 
%    % S(k,2).value = -out_maha{k}.solver.solver.objective;
%     S(k,1).value1 =  out_maha2.solver.solver.residual;
%   %  S(k,2).value1 = -out_maha2.solver.solver.objective;
%   
%     
%     
%     %S(2).value = out_osqp{k}.info.pri_res;
%     
%  
%     row = 2;
%     col = k;
%     for kk = 1:numel(S)
%         subplot(row,col,kk);
%         semilogy(S(kk).value);
%         hold on
%         semilogy(S(kk).value1);
%         hold on
%         semilogy(res_pri00);
%         hold on
%         semilogy(res_pri0);
%         hold on 
%         semilogy(res_pri10);
%         hold on 
%         semilogy(res_pri1)
%         hold on 
%         semilogy(res_pri20);
%         hold on
%         semilogy(res_pri2);
%         hold on
%         semilogy(res_pri30);
%         hold on
%         semilogy(res_pri3);
%         hold on 
%         legend('F-JNAMA','F-NAMA','NN-NAMA-QP','NA-NAMA-QP','JN-NAMA-QP','JA-NAMA-QP','RN-NAMA-QP','RA-NAMA-QP','DRN-NAMA-QP','DRA-NAMA-QP')
%         hold off
%         title(S(kk).site_name)
%     end
    
end