% NUMDIFF
% 
%   NUMDIFF(f, x, h), numerically computes the gradient of f at x

function J = numdiff(f, x, h)
    if nargin < 2, error('insufficient arguments'); end
    if nargin < 3, h = 1e-6; end
    n = length(x);
    v = f(x);
    m = length(v);
    J = zeros(m, n);
    for i=1:n
        x1 = x;
        x2 = x;
        x1(i) = x(i)+h;
        x2(i) = x(i)-h;
        v1 = f(x1);
        v2 = f(x2);
        J(:, i) = (v1-v2)'/(2*h);
    end
    % if f has values in R then transpose the Jacobian
    % (to get the gradient as a column-vector in this case)
    if m == 1
        J = J';
    end
end
