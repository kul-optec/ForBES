function obj = conjugate(f)
    if nargin < 1 || isempty(f)
        obj = [];
        return;
    end
    obj.isQuadratic = f.isConjQuadratic;
    obj.isConjQuadratic = f.isQuadratic;
    obj.hasHessian = 0;
    if isfield(f, 'makefconj')
        obj.makef = f.makefconj;
    end
    if isfield(f, 'makef')
        obj.makefconj = f.makef;
    end
    if isfield(f, 'makeprox')
        obj.makeprox = @() make_prox_conj(f);
    end
end

function fun = make_prox_conj(f)
    proxf = f.makeprox();
    fun = @(x, gam) call_prox_conj(x, gam, proxf);
end

function [p, v] = call_prox_conj(y, gam, prox)
    [z, v] = prox(y/gam, 1/gam);
    p = y - gam*z;
    v = y(:)'*z(:) - gam*(z(:)'*z(:)) - v;
end
