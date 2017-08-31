function [p, v] = prox(obj, x, gam)
    bx = obj.b .* x;
    ind = bx < 1;
    p(ind,1) = obj.b(ind) .* min(bx(ind)+gam*obj.mu,1);
    p(~ind,1) = x(~ind);
    v = obj.mu*sum(max(0,1-obj.b .* p));
end
