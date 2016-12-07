function [gam, cache, cache_bar, ops, flag] = Cache_BacktrackGamma(cache, gam, beta)

ops = Ops_Init();
flag = 0;
[gamma_ok, cache, cache_bar, ops1] = Cache_CheckGamma(cache, gam, beta);
ops = Ops_Sum(ops, ops1);
while ~gamma_ok
    gam = gam/2;
    flag = 1;
    [gamma_ok, cache, cache_bar, ops1] = Cache_CheckGamma(cache, gam, beta);
    ops = Ops_Sum(ops, ops1);
end

end