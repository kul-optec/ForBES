function [p, v] = prox(obj, x, gam)
% The implementation of the proximal mappings is based on:
% Cao, Sun, Xu, "Fast image deconvolution using closed-form
% thresholding formulas of L_q (q=1/2,2/3) regularization" (2013)
    mu = 2*gam*obj.lam;
    q = nthroot(54, 3)/4*nthroot(mu, 3)^2;
    absx = abs(x);
    phi = acos((mu/8)*sqrt((absx/3).^(-3)));
    ind0 = (absx <= q);
    p(ind0, 1) = 0;
    p(~ind0, 1) = (2/3)*(sign(x(~ind0)).*absx(~ind0).*(1+cos((2/3)*(pi-phi(~ind0)))));
    v = obj.lam*sum(abs(p).^(0.5));
end
