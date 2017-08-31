function [p, v] = prox(obj, x, gam)
    % use Moreau identity:
    %   prox(x; gam, f*) = x - gam*prox(x/gam; 1/gam, f)
    %
    [p1, v1] = obj.f.prox(x/gam, 1/gam);
    p = x - gam*p1;
    v = x(:)'*p1(:) - gam*p1(:)'*p1(:) - v1;
end
