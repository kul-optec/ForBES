function obj = diagOp(n, w)
    if nargin < 1, error('you should provide the dimension of the space'); end
    if numel(n) == 1, n = [n, 1]; end
    if nargin < 2, w = 1; end
    if numel(w) ~= 1 && any(size(w) ~= n), error('size of w must be 1 or n'); end
    obj.m = n;
    obj.n = n;
    obj.makeop = @() @(x) w.*x;
    obj.makeadj = @() @(y) w.*y;
end
