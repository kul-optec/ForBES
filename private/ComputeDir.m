function [dir, cacheDir] = ComputeDir(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir)
% in minfbe:    v = gradient of FBE
% in zerofpr:   v = residual

% compute direction according to the method
switch opt.methodID
    case {1, 7} % STEEPEST DESCENT and BARZILAI-BORWEIN
        [dir, cacheDir] = Direction_Steepest(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 2 % BFGS
        [dir, cacheDir] = Direction_BFGS(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 3 % L-BFGS
        [dir, cacheDir] = Direction_LBFGS(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 4 % CG-DESCENT
        [dir, cacheDir] = Direction_CG_descent(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 5 % CG-PRP
        [dir, cacheDir] = Direction_CG_PRP(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 6 % CG-DYHS
        [dir, cacheDir] = Direction_CG_DYHS(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 8 % Broyden
        [dir, cacheDir] = Direction_Broyden(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 9 % limited-memory Broyden (using compact representation)
        [dir, cacheDir] = Direction_LBroyden(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    case 10 % restarted Broyden (recursive formula)
        [dir, cacheDir] = Direction_RBroyden(prob, opt, it, hasGammaChanged, sk, yk, v, cacheDir);
    otherwise
        error('search direction not implemented');
end

dir = reshape(dir, prob.n);

end

function s = sgn(a)

s = 2*(a>=0)-1;

end
