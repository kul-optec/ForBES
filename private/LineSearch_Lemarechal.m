function [t, cachet, cachet1, ops, lsopt, exitflag] = LineSearch_Lemarechal(cache, dir, slope, t0, lsopt, varargin)
%LEMARECHALLS - computes a steplength t > 0 so that it satisfies the (weak) Wolfe conditions
%
% f(t) <= f(0) + delta*f'(0)
% f'(t) >= sigma*f'(0).
%
% exitflag = -1: gam is not small enough
% exitflag =  0: acceptable steplength was found
% exitflag =  1: maximum number of bracketing or sectioning iterations reached
% exitflag =  2: no further progress can be made
%
% Algorithm is described in Figure 1 of
%
% C. Lemarechal, A view of line searches, in: Optimization and Optimal Control, Auslender,
% Oettli, Stoer Eds, Lecture Notes in Control and Information Sciences 30,
% Springer Verlag (1981)
%
% see also
%
% J.-B. Hiriart-Urruty and C. Lemarechal (1996).
% Convex Analysis and Minimization Algorithms, vol I.
% Springer Verlag, Heidelberg, Algorithm 3.3.1 (Wolfe's line-search)

    % precompute stuff for the line search
    [cache, ops] = Cache_LineSearch(cache, dir);

    cachet1 = [];

    gam = cache.gam;

    t = t0;
    wolfe_hi = lsopt.delta*slope;
    wolfe_lo = lsopt.sigma*slope;
    a = 0; fa = cache.FBE; dfa = slope;
    tprev = a; fprev = fa; dfprev = dfa;
    b = inf; % upper bound
    rho = lsopt.rho;
    theta = lsopt.theta;
    exitflag = 1;
    testGammaFlag = 0;
    for it = 1:lsopt.nbracket
        [cachet, ops1] = Cache_LineFBE(cache, t, 1);
        ops = Ops_Sum(ops, ops1);
        if lsopt.testGamma && testGammaFlag
            [flagGamma, cachet, cachet1, ops1] = Cache_CheckGamma(cachet, gam, lsopt.beta);
            ops = Ops_Sum(ops, ops1);
            exitflag = flagGamma-1; % because CheckGamma returns 1 (good gamma) or 0 (bad gamma)
        end
        testGammaFlag = 0;
        if cachet.FBE > cache.FBE + t*wolfe_hi
            b = t; fb = cachet.FBE;
            if lsopt.interp == 1
                tn = LemarechalQuadInterp(a,fa,dfa,b,fb);
                % safeguard
                tn = min(tn,b - theta*(b - a));
                tn = max(tn,a + theta*(b - a));
            elseif lsopt.interp == 2
                [cachet, ops1] = Cache_LineFBE(cache, t, 3, cachet);
                ops = Ops_Sum(ops, ops1);
                dfb = cachet.dFBE;
                tn = LemarechalCubInterp(a,fa,dfa,b,fb,dfb);
                % safeguard
                tn = min(tn,b - theta*(b - a));
                tn = max(tn,a + theta*(b - a));
            else
                tn = 0.5*(a + b);
            end
            t = tn;
        else
            [cachet, ops1] = Cache_LineFBE(cache, t, 2, cachet);
            ops = Ops_Sum(ops, ops1);
            if cachet.dFBE < wolfe_lo
                a = t; fa = cachet.FBE; dfa = cachet.dFBE;
                if b == inf
                    % extrapolate
                    if lsopt.interp % we always have dfprev
                        tn = LemarechalCubInterp(tprev,fprev,dfprev,a,fa,dfa);
                        % safeguard
                        tn = max(tn,rho*tprev);
                    else
                        tn = rho*t;
                    end
                    testGammaFlag = 1;
                else
                    % interpolate
                    if lsopt.interp == 1
                        tn = LemarechalQuadInterp(a,fa,dfa,b,fb);
                        % safeguard
                        tn = min(tn,b - theta*(b - a));
                        tn = max(tn,a + theta*(b - a));
                    elseif lsopt.interp == 2
                        tn = LemarechalCubInterp(a,fa,dfa,b,fb,dfb);
                        % safeguard
                        tn = min(tn,b - theta*(b - a));
                        tn = max(tn,a + theta*(b - a));
                    else
                        tn = 0.5*(a + b);
                    end
                end
                tprev = t;fprev = fa;dfprev = dfa;
                t = tn;
            else
                exitflag = 0;
                break;
            end
        end

        if (b-a) <= lsopt.progTol
            exitflag = 2;
            break;
        end
    end

    if exitflag == 0 && lsopt.testGamma
        [flagGamma, cachet, cachet1, ops1] = Cache_CheckGamma(cachet, gam, lsopt.beta);
        ops = Ops_Sum(ops, ops1);
        exitflag = flagGamma-1; % because CheckGamma returns 1 (good gamma) or 0 (bad gamma)
    end
end

function t = LemarechalQuadInterp(t0,f0,df0,t1,f1)
    % Minimizer of interpolant belongs to [0,t1]
    q = f1-f0-t1*df0;
    q = 2*(f1-f0-(t1-t0)*df0)/(t1-t0)^2;
    if q > 0%quadratic is strongly convex
       c2 = df0-t0*q;
       t = -c2/q;
    else
        t = -1;
    end
end

function t = LemarechalCubInterp(t1,f1,df1,t2,f2,df2)
    % t1, t2 might not be sorted
    delta = t2 - t1 ;
    if  delta == 0
        t = t1;
    else
        d1 = df1 + df2 - 3*(f2-f1)/delta;
        d2 = d1^2 - df1*df2 ;
        if  d2 < 0 % /* complex roots, use secant method */
            if ( abs(df1) < abs (df2) )
                t = t1 - (t1 - t2)*(df1/(df1-df2)) ;
            elseif ( df1 ~= df1 )
                t = t2 - (t1 - t2)*(df2/(df1-df2)) ;
            else
                t = - 1 ;
            end
        else
            % first way: from Hager-Zhang code
    %                 d2 = sqrt(d2)*sign(delta);
    %                 v1 = df1 + d1 - d2 ;
    %                 v2 = df2 + d1 + d2 ;
    %                 if ( (v1 == 0) && (v2 == 0) )
    %                     t = -1;
    %                 elseif ( abs (v1) >= abs (v2) )
    %                     t = t1 + delta*df1/v1 ;
    %                 else
    %                     t = t2 - delta*df2/v2 ;
    %                 end
            %
            % second way: from Bonnans, Lemarechal
            d2 = sqrt(d2)*sign(delta);
            v1 = df2 + d2 - d1;
            v2 = df2 - df1 + 2*d2;
            if ( (v1 == 0) && (v2 == 0) )
                t = -1;
            else
                t = t2 - delta*(v1/v2);
            end
        end
    end
end
