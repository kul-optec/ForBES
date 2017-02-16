function [flag, cache_z] = Check_Gamma(cache, bet)

z = cache.Get_ProxGradStep();
cache_z = FBCache(cache.prob, z, cache.gam, cache.ops);
fz = cache_z.Get_f();
% fx = cache.Get_f();
% gradfx = cache.Get_Gradf();
% fprx = cache.Get_FPR();
if fz <= cache.fx - cache.gradfx(:)'*cache.FPR(:) + (1-bet)/(2*cache.gam)*(cache.normFPR()^2) + 1e-6*abs(cache.fx)
    flag = true;
else
    flag = false;
end
