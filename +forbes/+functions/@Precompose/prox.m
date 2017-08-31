function [p, v] = prox(obj, x, gam)
% See Prop. 23.32 in Bauschke, Combettes
% "Convex Analysis and Monotone Operator Theory in Hilber Spaces",
% 1st ed, 2011
    Lx = obj.L*x;
    [p1, v] = obj.f.prox(Lx, gam);
    p = x + obj.L'*(p1 - Lx)/obj.mu;
end
