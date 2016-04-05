function obj = hankelOp(p, q)
    if nargin < 2
        error('you must provide at least arguments p, q');
    end
    H = [];
    n = p+q-1;
    for i = 1:q
        H = [H; spdiags(ones(p,1), i-1, p, n)];
    end
    obj.m = p*q;
    obj.n = (p+q-1);
    obj.makeop = @() @(x) H*x;
    obj.makeadj = @() @(y) H'*y;
end
