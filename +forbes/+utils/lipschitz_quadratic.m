function Lf = lipschitz_quadratic(f, A, x)
    sizex = size(x);
    gradf_A0 = f.gradient(A*zeros(sizex));
    eigsOpt.issym = 1;
    eigsOpt.tol = 1e-3;
    funHessian = @(x) vec(A'*(f.gradient(A*reshape(x, sizex))-gradf_A0));
    Lf = eigs(funHessian, prod(sizex), 1, 'LM', eigsOpt);
end
