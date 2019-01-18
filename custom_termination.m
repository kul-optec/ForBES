function stop = custom_termination(prob, it, cache_0, cache_x,scale,Hessian)

 

 a = max(norm(-scale*cache_x.Get_Gradf(),'inf'), norm(scale*(cache_x.Get_Gradf() + cache_x.Get_FPR()/cache_x.Get_Gamma()), 'inf'));
% b = max(max(norm(Hessian*(cache_x.Get_Point()),'inf'),norm(prob.C1*cache_x.y,'inf')), norm(prob.q,'inf')); 
 rel_tol = 1e-3;
 abs_tol = 1e-3;
tol_primal = abs_tol + a*rel_tol;
%tol_dual = abs_tol + b*rel_tol;

%dual_residual = Hessian*(cache_x.Get_Point())+prob.C1*cache_x.y+prob.q;
normInfFPR = norm(cache_x.Get_FPR(), 'inf')/cache_x.Get_Gamma();
if normInfFPR <= tol_primal %& dual_residual <= tol_dual
    stop = true;
else
    stop = false;
end

end