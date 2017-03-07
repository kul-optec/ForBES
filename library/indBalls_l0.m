%INDBALLS_L0 Indicator function of the product of (nonconvex) L0 balls
% with given radius.
%
%   INDBALLS_L0(m, N)

function obj = indBalls_l0(m, N)
    obj.makeprox = @() @(x, gam) call_indBalls_l0_prox(x, m, N);
end

function [prox, val] = call_indBalls_l0_prox(x, m, N)
    n = length(x);
    k = n/m; % we are assuming n is a multiple of m
    prox = x;
    for i = 1:k
        [~, I] = sort(prox((i-1)*m+1:i*m), 'descend');
        prox((i-1)*m + I(N+1:end)) = 0;
    end
    val = 0;
end
