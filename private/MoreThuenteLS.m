% Copyright (C) 2015, Lorenzo Stella and Panagiotis Patrinos
%
% This file is part of ForBES.
% 
% ForBES is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% ForBES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with ForBES. If not, see <http://www.gnu.org/licenses/>.

function [cachet, t, cnt, info ] = MoreThuenteLS(prob, gam, cache, df0, lsopt)
%     Function MoreThuenteLS
%
%     The purpose of MoreThuenteLS is to find a step which satisfies
%     a sufficient decrease condition and a curvature condition.
%
%     At each stage MoreThuenteLS updates an interval of
%     uncertainty with endpoints stx and sty. The interval of
%     uncertainty is initially chosen so that it contains a
%     minimizer of the modified function
%
%          f(x+stp*s) - f(x) - delta*stp*(gradf(x)'s).
%
%     If a step is obtained for which the modified function
%     has a nonpositive function value and nonnegative derivative,
%     then the interval of uncertainty is chosen so that it
%     contains a minimizer of f(x+stp*s).
%
%     The algorithm is designed to find a step which satisfies
%     the sufficient decrease condition
%
%           f(x+stp*s) <= f(x) + delta*stp*(gradf(x)'s),
%
%     and the curvature condition
%
%           abs(gradf(x+stp*s)'s)) <= sigma*abs(gradf(x)'s).
%
%     If delta is less than sigma and if, for example, the function
%     is bounded below, then there is always a step which satisfies
%     both conditions. If no step can be found which satisfies both
%     conditions, then the algorithm usually stops when rounding
%     errors prevent further progress. In this case stp only
%     satisfies the sufficient decrease condition.
%
%
%	stp is a nonnegative variable. On input stp contains an
%         initial estimate of a satisfactory step. On output
%         stp contains the final estimate.
%
%   delta and sigma are nonnegative input variables. Termination
%         occurs when the sufficient decrease condition and the
%         directional derivative condition are satisfied.
%
%	progTol is a nonnegative input variable. Termination occurs
%         when the relative width of the interval of uncertainty
%	      is at most progTol.
%
%  tmin and tmax are nonnegative input variables which
%	  specify lower and upper bounds for the step.
%
%	maxfev is a positive integer input variable. Termination
%         occurs when the number of calls to fcn is at least
%         maxfev by the end of an iteration.
%
%	info is an integer output variable set as follows:
%
%	  info = 0  The sufficient decrease condition and the
%                   directional derivative condition hold.
%
%	  info = -2  Relative width of the interval of uncertainty
%		    is at most progTol.
%
%	  info = -1  Number of calls to fcn has reached maxfev.
%
%	  info = 4  The step is at the lower bound tmin.
%
%	  info = 5  The step is at the upper bound tmax.
%
%	  info = 6  Rounding errors prevent further progress.
%                   There may not be a step which satisfies the
%                   sufficient decrease and curvature conditions.
%                   Tolerances may be too small.
%
%       nf is an integer output variable set to the number of
%         calls to fcn.
%
%	wa is a work array of length n.
%
%     Subprograms called
%
%	user-supplied......fcn
%
%	MINPACK-supplied...cstep
%
%	FORTRAN-supplied...abs,max,min
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
%     **********
    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    nf = 0;
    p5 = .5;
    p66 = .66;
    xtrapf = 4;
    info = 0;
    infoc = 1;

    %     Compute the initial gradient in the search direction
    %     and check that s is a descent direction.
    %

    %
    %     Initialize local variables.
    %
    t = lsopt.tau0;
    f = cache.FBE;
    brackt = 0;
    stage1 = 1;
    finit = f;
    dgtest = lsopt.delta*df0;
    width = lsopt.tmax - lsopt.tmin;
    width1 = 2*width;
    %
    %     The variables stx, fx, dgx contain the values of the step,
    %     function, and directional derivative at the best step.
    %     The variables sty, fy, dgy contain the value of the step,
    %     function, and derivative at the other endpoint of
    %     the interval of uncertainty.
    %     The variables stp, f, dg contain the values of the step,
    %     function, and derivative at the current step.
    %
    stx = 0;
    fx = finit;
    dgx = df0;
    sty = 0;
    fy = finit;
    dgy = df0;
    %
    %     Start of iteration.
    %
    while (1)
        %
        %        Set the minimum and maximum steps to correspond
        %        to the present interval of uncertainty.
        %
        if (brackt)
            stmin = min(stx,sty);
            stmax = max(stx,sty);
        else
            stmin = stx;
            stmax = t + xtrapf*(t - stx);
        end
        %
        %        Force the step to be within the bounds tmax and tmin.
        %
        t = max(t,lsopt.tmin);
        t = min(t,lsopt.tmax);
        %
        %        If an unusual termination is to occur then let
        %        t be the lowest point obtained so far.
        %
        if ((brackt && (t <= stmin || t >= stmax)) || nf >= lsopt.maxfev-1 || infoc == 0 || (brackt && stmax-stmin <= lsopt.progTol*stmax))
            t = stx;
        end
        %
        %        Evaluate the function and gradient at t
        %        and compute the directional derivative.
        %
        [cachet, cnt1] = DirFBE(prob, gam, t, cache, 3);
        cnt = cnt+cnt1;
        nf = nf+1;
        f = cachet.FBE; dg = cachet.dFBE;

        ftest1 = finit + t*dgtest;
        %
        %        Test for convergence.
        %
        if ((brackt && (t <= stmin || t >= stmax)) || infoc == 0)
            % Rounding errors prevent further progress
            info = 6;
            return
        end
        if (t == lsopt.tmax && f <= ftest1 && dg <= dgtest)
            % The step is at the upper bound tmax
            info = 5;
            return
        end
        if (t == lsopt.tmin && (f > ftest1 || dg >= dgtest))
            % The step is at the lower bound tmin
            info = 4;
            return
        end
        if (nf >= lsopt.maxfev)
            info = -1;
            return
        end
        if (brackt && stmax-stmin <= lsopt.progTol*stmax)
            info = -2;
            return
        end
        if (f <= ftest1 && abs(dg) <= lsopt.sigma*(-df0))
            info = 0;
            return
        end
        %
        %        In the first stage we seek a step for which the modified
        %        function has a nonpositive value and nonnegative derivative.
        %
        if (stage1 && f <= ftest1 && dg >= min(lsopt.delta,lsopt.sigma)*df0)
            stage1 = 0;
        end
        %
        %        A modified function is used to predict the step only if
        %        we have not obtained a step for which the modified
        %        function has a nonpositive function value and nonnegative
        %        derivative, and if a lower function value has been
        %        obtained but the decrease is not sufficient.
        %
        if (stage1 && f <= fx && f > ftest1)
            %
            %           Define the modified function and derivative values.
            %
            fm = f - t*dgtest;
            fxm = fx - stx*dgtest;
            fym = fy - sty*dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;
            %
            %           Call cstep to update the interval of uncertainty
            %           and to compute the new step.
            %
            [stx,fxm,dgxm,sty,fym,dgym,t,fm,dgm,brackt,infoc] = MoreThuenteCstep(stx,fxm,dgxm,sty,fym,dgym,t,fm,dgm,brackt,stmin,stmax);
            %
            %           Reset the function and gradient values for f.
            %
            fx = fxm + stx*dgtest;
            fy = fym + sty*dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
        else
            %
            %           Call cstep to update the interval of uncertainty
            %           and to compute the new step.
            %
            [stx,fx,dgx,sty,fy,dgy,t,f,dg,brackt,infoc] = MoreThuenteCstep(stx,fx,dgx,sty,fy,dgy,t,f,dg,brackt,stmin,stmax);
        end
        %
        %        Force a sufficient decrease in the size of the
        %        interval of uncertainty.
        %
        if (brackt)
            if (abs(sty-stx) >= p66*width1)
                t = stx + p5*(sty - stx);
            end
            width1 = width;
            width = abs(sty-stx);
        end
        %
        %        End of iteration.
        %
    end
end

function  [stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,info] = MoreThuenteCstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)
%     Subroutine cstep
%
%     The purpose of cstep is to compute a safeguarded step for
%     a linesearch and to update an interval of uncertainty for
%     a minimizer of the function.
%
%     The parameter stx contains the step with the least function
%     value. The parameter stp contains the current step. It is
%     assumed that the derivative at stx is negative in the
%     direction of the step. If brackt is set true then a
%     minimizer has been bracketed in an interval of uncertainty
%     with endpoints stx and sty.
%
%     The subroutine statement is
%
%       subroutine cstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
%                        stpmin,stpmax,info)
%
%     where
%
%       stx, fx, and dx are variables which specify the step,
%         the function, and the derivative at the best step obtained
%         so far. The derivative must be negative in the direction
%         of the step, that is, dx and stp-stx must have opposite
%         signs. On output these parameters are updated appropriately.
%
%       sty, fy, and dy are variables which specify the step,
%         the function, and the derivative at the other endpoint of
%         the interval of uncertainty. On output these parameters are
%         updated appropriately.
%
%       stp, fp, and dp are variables which specify the step,
%         the function, and the derivative at the current step.
%         If brackt is set true then on input stp must be
%         between stx and sty. On output stp is set to the new step.
%
%       brackt is a logical variable which specifies if a minimizer
%         has been bracketed. If the minimizer has not been bracketed
%         then on input brackt must be set false. If the minimizer
%         is bracketed then on output brackt is set true.
%
%       stpmin and stpmax are input variables which specify lower
%         and upper bounds for the step.
%
%       info is an integer output variable set as follows:
%         If info = 1,2,3,4,5, then the step has been computed
%         according to one of the five cases below. Otherwise
%         info = 0, and this indicates improper input parameters.
%
%     Subprograms called
%
%       FORTRAN-supplied ... abs,max,min,sqrt
%                        ... dble
%
%     Argonne National Laboratory. MINPACK Project. June 1983
%     Jorge J. More', David J. Thuente
%
%     **********
    p66 = 0.66;
    info = 0;
    %
    %     Check the input parameters for errors.
    %
    if ((brackt && (stp <= min(stx,sty) || stp >= max(stx,sty))) || dx*(stp-stx) >= 0.0 || stpmax < stpmin)
        return
    end
    %
    %     Determine if the derivatives have opposite sign.
    %
    sgnd = dp*(dx/abs(dx));
    %
    %     First case. A higher function value.
    %     The minimum is bracketed. If the cubic step is closer
    %     to stx than the quadratic step, the cubic step is taken,
    %     else the average of the cubic and quadratic steps is taken.
    %
    if (fp > fx)
        info = 1;
        bound = 1;
        theta = 3*(fx - fp)/(stp - stx) + dx + dp;
        s = norm([theta,dx,dp],inf);
        gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
        if (stp < stx)
            gamma = -gamma;
        end
        p = (gamma - dx) + theta;
        q = ((gamma - dx) + gamma) + dp;
        r = p/q;
        stpc = stx + r*(stp - stx);
        stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx);
        if (abs(stpc-stx) < abs(stpq-stx))
            stpf = stpc;
        else
            stpf = stpc + (stpq - stpc)/2;
        end
        brackt = 1;
        %
        %     Second case. A lower function value and derivatives of
        %     opposite sign. The minimum is bracketed. If the cubic
        %     step is closer to stx than the quadratic (secant) step,
        %     the cubic step is taken, else the quadratic step is taken.
        %
    elseif (sgnd < 0.0)
        info = 2;
        bound = 0;
        theta = 3*(fx - fp)/(stp - stx) + dx + dp;
        s = norm([theta,dx,dp],inf);
        gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
        if (stp > stx)
            gamma = -gamma;
        end
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dx;
        r = p/q;
        stpc = stp + r*(stx - stp);
        stpq = stp + (dp/(dp-dx))*(stx - stp);
        if (abs(stpc-stp) > abs(stpq-stp))
            stpf = stpc;
        else
            stpf = stpq;
        end
        brackt = 1;
        %
        %     Third case. A lower function value, derivatives of the
        %     same sign, and the magnitude of the derivative decreases.
        %     The cubic step is only used if the cubic tends to infinity
        %     in the direction of the step or if the minimum of the cubic
        %     is beyond stp. Otherwise the cubic step is defined to be
        %     either stpmin or stpmax. The quadratic (secant) step is also
        %     computed and if the minimum is bracketed then the the step
        %     closest to stx is taken, else the step farthest away is taken.
        %
    elseif (abs(dp) < abs(dx))
        info = 3;
        bound = 1;
        theta = 3*(fx - fp)/(stp - stx) + dx + dp;
        s = norm([theta,dx,dp],inf);
        %
        %        The case gamma = 0 only arises if the cubic does not tend
        %        to infinity in the direction of the step.
        %
        gamma = s*sqrt(max(0.,(theta/s)^2 - (dx/s)*(dp/s)));
        if (stp > stx)
            gamma = -gamma;
        end
        p = (gamma - dp) + theta;
        q = (gamma + (dx - dp)) + gamma;
        r = p/q;
        if (r < 0.0 && gamma ~= 0.0)
            stpc = stp + r*(stx - stp);
        elseif (stp > stx)
            stpc = stpmax;
        else
            stpc = stpmin;
        end
        stpq = stp + (dp/(dp-dx))*(stx - stp);
        if (brackt)
            if (abs(stp-stpc) < abs(stp-stpq))
                stpf = stpc;
            else
                stpf = stpq;
            end
        else
            if (abs(stp-stpc) > abs(stp-stpq))
                stpf = stpc;
            else
                stpf = stpq;
            end
        end
        %
        %     Fourth case. A lower function value, derivatives of the
        %     same sign, and the magnitude of the derivative does
        %     not decrease. If the minimum is not bracketed, the step
        %     is either stpmin or stpmax, else the cubic step is taken.
        %
    else
        info = 4;
        bound = 0;
        if (brackt)
            theta = 3*(fp - fy)/(sty - stp) + dy + dp;
            s = norm([theta,dy,dp],inf);
            gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s));
            if (stp > sty)
                gamma = -gamma;
            end
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p/q;
            stpc = stp + r*(sty - stp);
            stpf = stpc;
        elseif (stp > stx)
            stpf = stpmax;
        else
            stpf = stpmin;
        end
    end
    %
    %     Update the interval of uncertainty. This update does not
    %     depend on the new step or the case analysis above.
    %
    if (fp > fx)
        sty = stp;
        fy = fp;
        dy = dp;
    else
        if (sgnd < 0.0)
            sty = stx;
            fy = fx;
            dy = dx;
        end
        stx = stp;
        fx = fp;
        dx = dp;
    end
    %
    %     Compute the new step and safeguard it.
    %
    stpf = min(stpmax,stpf);
    stpf = max(stpmin,stpf);
    stp = stpf;
    if (brackt && bound)
        if (sty > stx)
            stp = min(stx+p66*(sty-stx),stp);
        else
            stp = max(stx+p66*(sty-stx),stp);
        end
    end
end

