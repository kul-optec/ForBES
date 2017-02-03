close all;
clear;

ASSERT_EPS = 1e-14;

% two simple operators

diag1 = randn(5,1);
op1 = diagOp(5, diag1); % R^5 -> R^5
mat1 = diag(diag1);

mat2 = randn(5, 3);
op2 = matOp(mat2); % R^3 -> R^5

opsSum = sumOp({op1, op2}); % R^8 -> R^5
callOpsSum = opsSum.makeop();

for i = 1:100
    x = randn(8, 1);
    y1 = callOpsSum(x);
    y2 = [mat1, mat2]*x;
    assert(norm(y1-y2, inf) <= ASSERT_EPS);
end

callOpsSumAdj = opsSum.makeadj();

for i = 1:100
    y = randn(5, 1);
    x1 = callOpsSumAdj(y);
    x2 = [mat1'*y; mat2'*y];
    assert(norm(x1-x2, inf) <= ASSERT_EPS);
end

% two less simple operators

diag1 = randn(2,4);
op1 = diagOp([2,4], diag1); % R^{2x4} -> R^{2x4}

diag2 = randn(2,4);
op2 = diagOp([2,4], diag2); % R^{2x4} -> R^{2x4}

opsSum = sumOp({op1, op2}); % R^16 -> R^{2x4}
assert(all(opsSum.n == [16, 1]));
assert(all(opsSum.m == [2, 4]));
callOpsSum = opsSum.makeop();

for i = 1:100
    x = randn(16, 1);
    y1 = callOpsSum(x);
    y2 = diag1.*reshape(x(1:8), 2, 4) + diag2.*reshape(x(9:16), 2, 4);
    assert(norm(y1-y2, inf) <= ASSERT_EPS);
end

callOpsSumAdj = opsSum.makeadj();

for i = 1:100
    y = randn(2, 4);
    x1 = callOpsSumAdj(y);
    x2 = [diag1.*y, diag2.*y];
    x2 = x2(:);
    assert(norm(x1-x2, inf) <= ASSERT_EPS);
end
