% NUMDIFF
% 
%   NUMDIFF(f, x, h), numerically computes the gradient of f at x

function G = numdiff(f, x, h)
    if nargin < 2, error('insufficient arguments'); end
    if nargin < 3, h = 1e-6; end
    [n, k] = size(x);
    v = f(x);
    G = zeros(n, k);
    for i=1:n
        for j=1:k
            x1 = x;
            x2 = x;
            x1(i, j) = x(i, j)+h;
            x2(i, j) = x(i, j)-h;
            v1 = f(x1);
            v2 = f(x2);
            G(i, j) = (v1-v2)'/(2*h);
        end
    end
end
