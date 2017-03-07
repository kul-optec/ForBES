%INDBALL_L0 Indicator function of the (nonconvex) L0 ball with given radius.
%
%   INDBALL_L0(N) builds the function
%
%       g(x) = 0    if nnz(x) <= N
%            = +inf otherwise
%
%   Argument N is required.

function obj = indBall_l0(N)
    obj.makeprox = @() @(x, gam) call_indBall_l0_prox(x, N);
end

function [prox, val] = call_indBall_l0_prox(x, N)
    prox = x;
    [~, I] = sort(prox, 'descend');
    prox(I(N+1:end)) = 0;
    val = 0;
end
