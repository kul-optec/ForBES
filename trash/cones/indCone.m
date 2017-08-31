% This is the function
%
%   f(x) = ind_K(x)
%
% where K is a cone in SCS format.

function obj = indCone(K)
    K = validate_cone(K);
    obj.isConvex = 1;
    obj.makeprox = @() @(z, gam) call_indCone_proj(z, K);
end

function [proj, val] = call_indCone_proj(z, K)
    proj = proj_cone(z, K);
    val = 0;
end
