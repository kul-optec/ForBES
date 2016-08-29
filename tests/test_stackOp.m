ASSERT_EPS = 1e-14;

% two simple operators

diag1 = randn(5,1);
op1 = diagOp(5, diag1); % R^5 -> R^5
mat1 = diag(diag1);

mat2 = randn(3, 5);
op2 = matOp(mat2); % R^5 -> R^3

opsStack = stackOp({op1, op2}); % R^5 -> R^8
callOpsStack = opsStack.makeop();

for i = 1:100
    x = randn(5, 1);
    y1 = callOpsStack(x);
    y2 = [mat1*x; mat2*x];
    assert(norm(y1-y2, inf) <= ASSERT_EPS);
end

callOpsStackAdj = opsStack.makeadj();

for i = 1:100
    y = randn(8, 1);
    x1 = callOpsStackAdj(y);
    x2 = [mat1', mat2']*y;
    assert(norm(x1-x2, inf) <= ASSERT_EPS);
end

% two less simple operators

diag1 = randn(2,4);
op1 = diagOp([2,4], diag1); % R^{2x4} -> R^{2x4}

diag2 = randn(2,4);
op2 = diagOp([2,4], diag2); % R^{2x4} -> R^{2x4}

opsStack = stackOp({op1, op2}); % R^{2x4} -> R^{16}
assert(all(opsStack.n == [2, 4]));
assert(all(opsStack.m == [16, 1]));
callOpsStack = opsStack.makeop();

for i = 1:100
    x = randn(2, 4);
    y1 = callOpsStack(x);
    y2 = [vec(diag1.*x); vec(diag2.*x)];
    assert(norm(y1-y2, inf) <= ASSERT_EPS);
end

callOpsStackAdj = opsStack.makeadj();

for i = 1:100
    y = randn(16, 1);
    x1 = callOpsStackAdj(y);
    x2 = reshape(diag1(:).*y(1:8) + diag2(:).*y(9:16), 2, 4);
    assert(norm(x1-x2, inf) <= ASSERT_EPS);
end
