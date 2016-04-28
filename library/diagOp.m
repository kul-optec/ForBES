function obj = diagOp(n, w)
    if nargin < 1
        error('you should provide the dimension of the space');
    end
    if nargin < 2
        w = 1;
    end
    if length(w) ~= 1 && length(w) ~= n
        error('length of w must be 1 or n');
    end
    obj.m = n;
    obj.n = n;
    obj.makeop = @() @(x) w.*x;
    obj.makeadj = @() @(y) w.*y;
end
