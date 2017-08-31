% This is the function
%
%   f(w) = (b/2)||(x,y)-a||^2 + ind_K*(y)
%
% where K* is the dual cone to K and x has length n.
% Only the conjugate function is implemented here.

function obj = smoothCP(n, K, a, b)
    K = validate_cone(K);
    obj.isConvex = 1;
    obj.makefconj = @() @(z) call_coneCost_conj(z, n, K, a, b);
end

function [v, g] = call_coneCost_conj(z, n, K, a, b)
    p1 = a(1:n) + z(1:n)/b;
    p2 = proj_dual_cone(a(n+1:end) + z(n+1:end)/b, K);
    x = z(1:n);
    y = proj_dual_cone(z(n+1:end), K);
    g = [p1; p2];
    v = z'*g - (b/2)*norm(g-a)^2;
end
