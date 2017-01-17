% This is the function
%
%   f(x) = ind_S(x)
%
% where S = 0^n \times K \times 0, and K is a cone in SCS format.

function obj = nonsmoothCP(m, n, K)
    K = validate_cone(K);
    obj.isConvex = 1;
    obj.makeprox = @() @(z, gam) call_indCone_proj(z, m, n, K);
end

function [proj, val] = call_indCone_proj(z, m, n, K)
    proj = [zeros(n, 1); proj_cone(z(n+1:n+m), K); 0];
    val = 0;
end
