%L2NORM Allocates the L2 norm function.
%
%   L2NORM(mu) builds the function
%
%       g(x) = mu*||x||_2

function obj = l2Norm(mu)
    if nargin < 1
        mu = 1;
    end
    obj.makeprox = @() @(x, gam) call_l2Norm_prox(x, gam, mu);
end

function [prox, val] = call_l2Norm_prox(x, gam, mu)
    normx = sqrt(x'*x);
    mugam = mu*gam;
    if normx <= mugam
        prox = zeros(length(x),1);
        val = 0;
    else
        scal = (1-mugam/normx);
        prox = (1-mugam/normx)*x;
        val = mu*scal*normx;
    end
end
