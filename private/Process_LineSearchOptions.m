function lsopt = Process_LineSearchOptions(opt)

%     %  factor in [0, 1] used to compute average cost magnitude C_k as follows:
%     % Q_k = 1 + (Delta)Q_k-1, Q_0 = 0,  C_k = C_k-1 + (|f_k| - C_k-1)/Q_k
%     lsopt.Delta = 0.7;% this goes here to include Hager-Zhang line search as a backup
    % Armijo condition parameter delta, range [0, .5]
    % phi (a) - phi (0) <= delta phi'(0)
    lsopt.testGamma = opt.adaptive && ~strcmp(opt.solver, 'zerofpr');
    if isfield(opt, 'delta'), lsopt.delta = opt.delta;
    else lsopt.delta = 0.1; end
    lsopt.beta = opt.beta;
    switch opt.linesearch
        case ''
            lsopt.linesearchfun = @() error('no line-search selected');
        case 'backtracking'
            lsopt.linesearchfun = @LineSearch_Backtracking;
            lsopt.progTol = 0;
            lsopt.nLS = 50;
        case 'backtracking-nm' % nonmonotone
            lsopt.linesearchfun = @LineSearch_BacktrackingNM;
            lsopt.eta = 0.85;
            lsopt.progTol = 0;
            lsopt.nLS = 50;
        case 'backtracking-armijo'
            lsopt.linesearchfun = @LineSearch_BacktrackingArmijo;
            lsopt.progTol = 0;
            lsopt.nLS = 50;
        case 'lemarechal'
            lsopt.linesearchfun = @LineSearch_Lemarechal;
            lsopt.sigma = 0.9;
            % maximum number of iterations
            lsopt.nbracket = 100;
            % type of interpolation - 0 [bisection], 1 [quadratic
            % interpolation], 2 [cubic interpolation when possible]
            lsopt.interp = 0;
            if isfield(opt, 'interp'), lsopt.interp = opt.interp; end
            % stop when length of interval is below progTol
            lsopt.progTol = 0;
            %  growth factor in search for initial bracket interval
            lsopt.rho = 5;
            % parameter for safe-guarding (must be in (0,1/2])
            lsopt.theta = 0.49;
        case 'hager-zhang'
            lsopt.linesearchfun = @LineSearch_HagerZhang;
            lsopt.sigma = 0.9;
            %  maximum number of times the bracketing interval grows during expansion
            lsopt.nexpand = 50;
            % maximum number of secant steps
            lsopt.nsecant = 50;
            % maximum number of times the bracketing interval contracts
            lsopt.ncontract = 10;
            % factor by which eps grows when line search fails during contraction
            lsopt.egrow = 10;
            lsopt.QuadOK = true;
            % T => when possible, use a cubic step in the line search
            lsopt.UseCubic = true;
            % true => estimated error in function value is eps*Ck,
            % false => estimated error in function value is eps */
            lsopt.PertRule = true;
            lsopt.eps = 1e-6;
            % |f| < SmallCost*starting cost => skip QuadStep and set PertRule = FALSE*/
            lsopt.SmallCost = 1e-30;
            % T => use approximate Wolfe line search
            % F => use ordinary Wolfe line search, switch to approximate Wolfe when
            % |f_k+1-f_k| < omega*C_k, C_k = average size of cost */
            lsopt.AWolfe = false;
            lsopt.omega = 1e-3;
            % factor by which secant step is amplified during expansion phase where minimizer is bracketed
            lsopt.SecantAmp = 1.05;
            % factor by which rho grows during expansion phase where minimizer is bracketed
            lsopt.RhoGrow = 2.0;
            %  maximum number of times that eps is updated
            lsopt.neps = 5;
            % maximum factor secant step increases stepsize in expansion phase
            lsopt.ExpandSafe = 200;
            % value of the parameter theta in the cg_descent update formula:
            % W. W. Hager and H. Zhang, A survey of nonlinear conjugate gradient
            % methods, Pacific Journal of Optimization, 2 (2006), pp. 35-58.
            lsopt.theta = 0.5;
            %  growth factor in search for initial bracket interval
            lsopt.rho = 5;
            % decay factor for bracket interval width in line search, range (0, 1)
            lsopt.gamma = 0.66;
        case 'more-thuente'
            lsopt.linesearchfun = @LineSearch_MoreThuente;
            lsopt.sigma = 0.9;
            lsopt.progTol = 0;
            lsopt.tmin = 0;
            lsopt.tmax = 1e15;
            lsopt.maxfev = 100;
        case 'fletcher'
            lsopt.linesearchfun = @LineSearch_Fletcher;
            lsopt.sigma = 0.9;
            %  maximum number of times the bracketing interval grows during expansion
            lsopt.nbracket = 50;
            % maximum number of section steps
            lsopt.nsection = 50;
            % stop when progress is below progTol
            lsopt.progTol = 0;
            % estimate of minimum value of the function
            lsopt.fmin = -inf;
        otherwise
            error('unknown line-search');
    end

    % if method is not L-BFGS then initial stepsize is selected according to Hager-Zhang
    lsopt.quadStep = true;
    % starting guess for line search =
    % psi0 ||x_0||_infty over ||g_0||_infty if x_0 != 0
    % psi0 |f(x_0)|/||g_0||_2               otherwise */
    lsopt.psi0 = 0.01;
    % when the function is approximately quadratic, use gradient at
    % psi1*psi2*previous step for estimating initial stepsize */
    lsopt.psi1 = 1.0 ;
    % when starting a new cg iteration, our initial guess for the line
    % search stepsize is psi2*previous step */
    lsopt.psi2 = 2;

end
