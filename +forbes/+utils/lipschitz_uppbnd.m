function Lf = lipschitz_uppbnd(f, A, x)
    delta = max(1e-12, x*1e-6);
    y = x + delta;
    [gradf_Ax, ~] = f.gradient(A*x);
    [gradf_Ay, ~] = f.gradient(A*y);
    Lf = norm(A'*(gradf_Ax - gradf_Ay), 'fro')/norm(delta, 'fro');
end
