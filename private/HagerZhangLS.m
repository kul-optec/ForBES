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

function [cachet, alpha, cnt, info] = HagerZhangLS(prob, gam, cache, df0, lsopt)
% Hager Zhang line search based on CG-DESCENT Version 6.7  (April 7, 2014)
%    Approximate Wolfe line search routine
%    info:
%        0 (Wolfe or approximate Wolfe conditions satisfied)
%        3 (slope always negative in line search)
%        4 (number of line search iterations exceed nline)
%        6 (excessive updating of eps)
%        7 (Wolfe conditions never satisfied)
%    ========================================================================= */

    AWolfe = lsopt.AWolfe ;
    alpha = lsopt.tau0;
    f0 = cache.FBE;
    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    
    if lsopt.PertRule
        epsilon = lsopt.eps*abs(f0);
    else
        epsilon = lsopt.eps;
    end
    
    lsopt.wolfe_hi  = lsopt.delta*df0;
    lsopt.wolfe_lo  = lsopt.sigma*df0;
    lsopt.awolfe_hi = (2*lsopt.delta - 1)*df0;
    lsopt.fpert     = f0 + epsilon;

    % evaluate function or gradient at alpha (starting guess)
    if ( lsopt.QuadOK )
        [cachet, cnt1] = DirFBE(prob, gam, alpha, cache, 3);
        cnt = cnt+cnt1;
        f = cachet.FBE; df = cachet.dFBE;
        fb = f;
        if ( ~AWolfe ), fb = fb - alpha*lsopt.wolfe_hi ;end
        qb = true ; % function value at b known
    else
        [cachet, cnt1] = DirFBE(prob, gam, alpha, cache, 2);
        cnt = cnt+cnt1;
        df = cachet.dFBE;
        qb = false ;
    end
    b = alpha ;

    if ( AWolfe )
        db = df ;
        d0 = df0 ;
        da = df0;
    else
        db = df - lsopt.wolfe_hi ;
        d0 = df0 - lsopt.wolfe_hi ;
        da = d0;
    end
    a  = 0 ;
    a1 = 0 ;
    d1 = d0 ;
    fa = f0 ;

    %     if a quadratic interpolation step performed, check Wolfe conditions */
    if ( (lsopt.QuadOK) && (f <= f0) )
        if ( HagerZhangTestWolfe (alpha, f, df, f0,lsopt) ), info = 0; return ;end
    end

    %       Find initial interval [a,b] such that
    %       da <= 0, db >= 0, fa <= fpert = [(f0 + eps*abs(f0)) or (f0 + eps)] */
    rho = lsopt.rho ;
    ngrow = 1 ;
    while ( db < 0 )
        if ( ~qb )
            [cachet, cnt1] = DirFBE(prob, gam, alpha, cache, 1, cachet);
            cnt = cnt+cnt1;
            f = cachet.FBE;
            if ( AWolfe )
                fb = f ;
            else
                fb = f - b*lsopt.wolfe_hi ;
            end
            qb = true ;
        end
        if ( fb > lsopt.fpert ) % contract interval [a, b]
            [a,fa,da,b,fb,db,alpha,status,lsopt,cachet,cnt1] = HagerZhangUpdate (a, fa, da, b, fb, db,prob,gam,cache,lsopt,f0) ;
            cnt = cnt+cnt1;
            if ( status == 0 ), info = 0; return ;end  % /* Wolfe conditions hold */
            if ( status == -2 ), break,end ; %/* db >= 0 */
            if ( lsopt.neps > 0 ), info = 6;return,end
        end

        % expansion phase 
        ngrow = ngrow +1 ;
        if ( ngrow > lsopt.nexpand ), info = 3; return,end
        % update interval (a replaced by b) */
        a = b ;
        fa = fb ;
        da = db ;
        % store old values of a and corresponding derivative
        d2 = d1 ;
        d1 = da ;
        a2 = a1 ;
        a1 = a ;

        bmin = rho*b ;
        if ( (ngrow == 2) || (ngrow == 3) || (ngrow == 6) )
            if ( d1 > d2 )
                if ( ngrow == 2 )
                    b = a1 - (a1-a2)*(d1/(d1-d2)) ;
                else
                    if ( (d1-d2)/(a1-a2) >= (d2-d0)/a2 )
                        % convex derivative, secant overestimates minimizer
                        b = a1 - (a1-a2)*(d1/(d1-d2)) ;
                    else
                        % concave derivative, secant underestimates minimizer
                        b = a1 - lsopt.SecantAmp*(a1-a2)*(d1/(d1-d2)) ;
                    end
                end
                % safeguard growth 
                b = min (b, lsopt.ExpandSafe*a1) ;
            else
                rho = rho*lsopt.RhoGrow ;
            end
        else
            rho = rho*lsopt.RhoGrow ;
        end
        b = max (bmin, b) ;
        alpha = b ;
        [cachet, cnt1] = DirFBE(prob, gam, alpha, cache, 2);
        cnt = cnt+cnt1;
        df = cachet.dFBE;
        b = alpha ;
        qb = false ;
        if ( AWolfe )
            db = df ;
        else
            db = df - lsopt.wolfe_hi ;
        end

    end

    %     /* we now have fa <= fpert, da >= 0, db <= 0 */
    toggle = 0 ;
    width = b - a ;
    qb0 = false ;
    for iter = 0:lsopt.nsecant-1
        %         /* determine the next iterate */
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) )
            lsopt.QuadOK = true ;
            if ( lsopt.UseCubic && qb )
                alpha = HagerZhangCubInterp (a, fa, da, b, fb, db) ;
                if ( alpha < 0 ) %/* use secant method */
                    if ( -da < db )
                        alpha = a - (a-b)*(da/(da-db)) ;
                    elseif ( da ~= db )
                        alpha = b - (a-b)*(db/(da-db)) ;
                    else
                        alpha = -1. ;
                    end
                end
            else
                if  ( -da < db )
                    alpha = a - (a-b)*(da/(da-db)) ;
                elseif ( da ~= db )
                    alpha = b - (a-b)*(db/(da-db)) ;
                else
                    alpha = -1. ;
                end
            end
            width = lsopt.gamma*(b - a) ;

        elseif ( toggle == 1 ) %/* iteration based on smallest value*/
            lsopt.QuadOK = true ;
            if ( lsopt.UseCubic )
                if ( alpha == a ) %/* a is most recent iterate */
                    alpha = HagerZhangCubInterp (a0, fa0, da0, a, fa, da) ;
                elseif ( qb0 )        %/* b is most recent iterate */
                    alpha = HagerZhangCubInterp (b, fb, db, b0, fb0, db0) ;
                else
                    alpha = -1. ;
                end
                %                 /* if alpha no good, use cubic between a and b */
                if ( (alpha <= a) || (alpha >= b) )
                    if ( qb )
                        alpha = HagerZhangCubInterp (a, fa, da, b, fb, db) ;
                    else
                        alpha = -1. ;
                    end
                end

                %                 /* if alpha still no good, use secant method */
                if ( alpha < 0 )
                    if ( -da < db )
                        alpha = a - (a-b)*(da/(da-db)) ;
                    elseif ( da ~= db )
                        alpha = b - (a-b)*(db/(da-db)) ;
                    else
                        alpha = -1. ;
                    end
                end
            else %/* ( use secant ) */
                if ( (alpha == a) && (da > da0) ) %/* use a0 if possible */
                    alpha = a - (a-a0)*(da/(da-da0)) ;
                elseif ( db < db0 )                   %/* use b0 if possible */
                    alpha = b - (b-b0)*(db/(db-db0)) ;
                else %/* secant based on a and b */
                    if      ( -da < db )
                        alpha = a - (a-b)*(da/(da-db)) ;
                    elseif ( da ~= db )
                        alpha = b - (a-b)*(db/(da-db)) ;
                    else
                        alpha = -1. ;
                    end
                end            
                if ( (alpha <= a) || (alpha >= b) )
                    if      ( -da < db )
                        alpha = a - (a-b)*(da/(da-db)) ;
                    elseif ( da ~= db )
                        alpha = b - (a-b)*(db/(da-db)) ;
                    else
                        alpha = -1. ;
                    end
                end
            end
        else
            alpha = .5*(a+b) ; %/* use bisection if b-a decays slowly */
            lsopt.QuadOK = false ;
        end
        if ( (alpha <= a) || (alpha >= b) )
            alpha = .5*(a+b) ;
            if ( (alpha == a) || (alpha == b) ),
                info = 7; 
                return ;
            end
            lsopt.QuadOK = false ; %/* bisection was used */
        end

        if ( toggle == 0 ) %/* save values for next iteration */
            a0 = a ;
            b0 = b ;
            da0 = da ;
            db0 = db ;
            fa0 = fa ;
            if ( qb )
                fb0 = fb ;
                qb0 = true ;
            end
        end

        toggle = toggle + 1 ;
        if ( toggle > 2 ), toggle = 0 ;end
        [cachet, cnt1] = DirFBE(prob, gam, alpha, cache, 3);
        cnt = cnt+cnt1;
        f = cachet.FBE;df = cachet.dFBE;
        if ( lsopt.QuadOK )
            if ( HagerZhangTestWolfe (alpha, f, df, f0,lsopt) )
                info = 0;
                return
            end
        end

        if ( ~AWolfe )
            f = f - alpha*lsopt.wolfe_hi ;
            df = df - lsopt.wolfe_hi ;
        end

        if ( df >= 0 )
            b = alpha ;
            fb = f ;
            db = df ;
            qb = true ;
        elseif ( f <= lsopt.fpert )
            a = alpha ;
            da = df ;
            fa = f ;
        else
            B = b ;
            if ( qb ), fB = fb ;end
            dB = db ;
            b = alpha ;
            fb = f ;
            db = df ;
            %             /* contract interval [a, alpha] */
            [a,fa,da,b,fb,db,alpha,status,lsopt,cachet,cnt1] = HagerZhangUpdate (a,fa,da,b,fb,db,prob,gam,cache,lsopt,f0) ;
            cnt = cnt+cnt1;
            if ( status == 0 ), info = 0; return; end
            if ( status == -1 ) %/* eps reduced, use [a, b] = [alpha, b] */
                if ( lsopt.neps > 5 ), info = 6; return; end
                a = b ;
                fa = fb ;
                da = db ;
                b = B ;
                if ( qb ), fb = fB ;end
                db = dB ;
            else
                qb = true ;
            end
        end
    end
    info = 4;
end

% /* =========================================================================
%    ==== update ========================================================
%    =========================================================================
%    The input for this routine is an interval [a, b] with the property that
%    fa <= fpert, da >= 0, db >= 0, and fb >= fpert. The returned status is
%
%   11  function or derivative not defined
%    0  if the Wolfe conditions are satisfied
%   -1  if a new value for eps is generated with the property that for the
%       corresponding fpert, we have fb <= fpert
%   -2  if a subinterval, also denoted [a, b], is generated with the property
%       that fa <= fpert, da >= 0, and db <= 0
%
%    NOTE: The input arguments are unchanged when status = -1
%    ========================================================================= */
function [a,fa,da,b,fb,db,alpha,info,lsopt,cachet,cnt] = HagerZhangUpdate(a,fa,da,b,fb,db,prob,gam,cache,lsopt,f0)
    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    AWolfe = lsopt.AWolfe ;
    f1 = fb ;
    toggle = 0 ;
    width = 0 ;
    for iter=0:lsopt.ncontract
        if ( (toggle == 0) || ((toggle == 2) && ((b-a) <= width)) )
            %             /* cubic based on bracketing interval */
            alpha = HagerZhangCubInterp (a, fa, da, b, fb, db) ;
            toggle = 0 ;
            width = lsopt.gamma*(b-a) ;
            if ( iter ), 
                lsopt.QuadOK = true ;
            end %/* at least 2 cubic iterations */
        elseif ( toggle == 1 )
            lsopt.QuadOK = true ;
            %             /* cubic based on most recent iterate and smallest value */
            if ( old < a ) %/* a is most recent iterate */
                alpha = HagerZhangCubInterp (a, fa, da, old, fold, dold) ;
            else         %  /* b is most recent iterate */
                alpha = HagerZhangCubInterp (a, fa, da, b, fb, db) ;
            end
        else
            alpha = .5*(a+b) ; %/* use bisection if b-a decays slowly */
            lsopt.QuadOK = false ;
        end

        if ( (alpha <= a) || (alpha >= b) )
            alpha = .5*(a+b) ;
            lsopt.QuadOK = false ; %/* bisection was used */
        end
        toggle = toggle + 1 ;
        if ( toggle > 2 ) 
            toggle = 0 ;
        end
        [cachet, cnt] = DirFBE(prob, gam, alpha, cache, 3);
        f = cachet.FBE; df = cachet.dFBE;

        if ( lsopt.QuadOK )
            if ( HagerZhangTestWolfe (alpha, f, df, f0, lsopt) ), 
                info = 0; 
                return 
            end 
        end

        if ( ~AWolfe )
            f = f - alpha*lsopt.wolfe_hi ;
            df = df - lsopt.wolfe_hi ;
        end
        if ( df >= 0 )
            a = alpha ;
            fb = f ;
            db = df ;
            info = -2;
            return
        end
        if ( f <= lsopt.fpert ) %/* update a using alpha */
            old = a ;
            fold = fa ;
            dold = da ;
            a = alpha ;
            fa = f ;
            da = df ;
        else                     %/* update b using alpha */
            old = b ;
            fold = fb;
            dold = db;
            b = alpha ;
            fb = f ;
            db = df ;
        end

    end

    %% This might need debugging
    %   see if the cost is small enough to change the PertRule 
    if ( abs (fb) <= lsopt.SmallCost ),
        lsopt.PertRule = false ;
    end

    %  increase eps if slope is negative after Parm->nshrink iterations 
    if ( lsopt.PertRule )
        if ( f0 ~= 0)
            lsopt.eps = lsopt.egrow*(f1-f0)/abs (f0) ;
            lsopt.fpert = f0 + abs (f0)*lsopt.eps ;
        else
            lsopt.fpert = 2*f1 ;
        end
    else
        lsopt.eps = lsopt.egrow*(f1-f0) ;
        lsopt.fpert = f0 + lsopt.eps ;
    end
    lsopt.neps = lsopt.neps+1 ;
    info = -1 ;
end

function done = HagerZhangTestWolfe (alpha, f, df,f0,lsopt)
    done = false;
    if ( df >= lsopt.wolfe_lo )
        % c test original Wolfe conditions
        if ( f-f0 <= alpha*lsopt.wolfe_hi )
            done = true;
            % c test approximate Wolfe conditions
        elseif ( lsopt.AWolfe )
            done = ( (f <= lsopt.fpert) & (df <= lsopt.awolfe_hi));
        end
    end
end

function t = HagerZhangCubInterp(t1,f1,df1,t2,f2,df2)
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
            d2 = sqrt(d2)*sign(delta);
            v1 = df1 + d1 - d2 ;
            v2 = df2 + d1 + d2 ;
            if ( (v1 == 0) && (v2 == 0) )
                t = -1;
            elseif ( abs (v1) >= abs (v2) )
                t = t1 + delta*df1/v1 ;
            else
                t = t2 - delta*df2/v2 ;
            end
            %
            % second way: from Bonnans, Lemarechal
            %         d2 = sqrt(d2)*sign(delta);
            %         v1 = df2 + d2 - d1;
            %         v2 = df2 - df1 + 2*d2;
            %         if ( (v1 == 0) && (v2 == 0) )
            %             t = -1;
            %         else
            %             t = t2 - delta*(v1/v2);
            %         end
        end
    end
end