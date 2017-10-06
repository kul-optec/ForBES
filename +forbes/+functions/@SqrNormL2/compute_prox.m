function [prox, val] = compute_prox(obj, x, gam)
    wgam = obj.w .* gam;
    prox = x./(1+wgam);
    val = ((obj.w .* prox(:))'*prox(:))/2;
end
