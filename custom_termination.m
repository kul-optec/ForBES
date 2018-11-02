function stop = custom_termination(prob, it, cache_0, cache_x,scale)

a = max(norm(-scale*cache_x.Get_Gradf(),'inf'), norm(scale*(cache_x.Get_Gradf() + cache_x.Get_FPR()/cache_x.Get_Gamma()), 'inf'));
%b = max(norm((cache_x.Get_Gradf()-prob.q),'inf'),norm(cache_x.Get_Gradf(),'inf'), norm(prob.q,'inf')); 
rel_tol = 1e-3;
abs_tol = 1e-3;
tol = abs_tol + a*rel_tol;



normInfFPR = norm(scale*cache_x.Get_FPR(), 'inf')/cache_x.Get_Gamma();
if normInfFPR <= tol
    stop = true;
else
    stop = false;
end

end