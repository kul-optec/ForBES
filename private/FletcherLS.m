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

function [cachet, t, cnt, exitflag ] = FletcherLS(prob,gam,cache,df0,lsopt)
%FletcherLS - computes a steplength t > 0 so that it satisfies the strong Wolfe conditions
%
% f(t) <= f(0) + delta*f'(0) 
% abs(f'(t)) <= -sigma*f'(0).
%
% exitflag =  0: acceptable steplength was found
% exitflag =  1: steplength t for which f(t) < fminimum was found 
% exitflag = -1: maximum number of bracketing or sectioning iterations reached
% exitflag = -2: no further progress can be made
%
% Algorithm is described in
%
% R. Fletcher, Practical Methods of Optimization, John Wiley & Sons, 1987,
% second edition, section 2.6.
   
    lsopt.wolfe_hi  = lsopt.delta*df0;
    lsopt.wolfe_lo  = lsopt.sigma*df0;
    f0 = cache.FBE;
    t0 = lsopt.tau0;
    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];

    % Find a bracket of acceptable points
    [cachet,a,b,fa,dfa,fb,dfb,t,exitflag,cnt1] = FletcherBracket(cache,prob,gam,f0,df0,t0,lsopt);
    cnt = cnt+cnt1;

    if exitflag == 2 
      % BracketingPhase found a bracket containing acceptable points; now find acceptable point 
      % within bracket
        [cachet,t,exitflag,cnt1] = FletcherSection(cache,prob,gam,f0,a,fa,dfa,b,fb,dfb,lsopt);
        cnt = cnt+cnt1;
    end
end

%-----------------------------------------------------------------------------------
function [cachet,a,b,fa,dfa,fb,dfb,t,exitflag,cnt] = FletcherBracket(cache,prob,gam,f0,df0,t0,lsopt) 
% 
% bracketingPhase finds a bracket [a,b] that contains acceptable points; a bracket 
% is the same as a closed interval, except that a > b is allowed.
%
% The outputs fa and dfa are the values of the function and the derivative 
% evaluated at the bracket endpoint 'a'. Similar notation applies to the endpoint 
% 'b'. The possible values of exitflag are like in LINESEARCH, with the additional 
% value exitflag = 2, which indicates that a bracket containing acceptable points 
% was found.

    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    tau1 = 9; % factor to expand the current bracket
    a = []; b = []; fa = []; dfa = []; fb = []; dfb = [];
    ft = f0; dft = df0;    

    % Set maximum value of t (determined by fminimum)
    tmax = (lsopt.fmin - f0)/(lsopt.wolfe_hi); 
    told = 0;

    % First trial t is user-supplied
    t = t0;
    for nbracket = 1:lsopt.nbracket 
        fold = ft; dfold = dft;
        [cachet, cnt1] = DirFBE(prob, gam, t, cache, 3);
        cnt = cnt+cnt1;
        ft = cachet.FBE; dft = cachet.dFBE;

        % Terminate if f < fminimum
        if ft <= lsopt.fmin
            exitflag = 1;
            return 
        end

        % Bracket located - case 1
        if ft > f0 + t*lsopt.wolfe_hi  || ft >= fold    
            a = told; fa = fold; dfa = dfold;
            b = t;    fb = ft;   dfb = dft;
            exitflag = 2;
            return 
        end

        % Acceptable steplength found; no need to call sectioning phase
        if abs(dft) <= -lsopt.wolfe_lo
            exitflag = 0;
            return
        end

        % Bracket located - case 2  
        if dft >= 0
            a = t;    fa = ft;   dfa = dft;
            b = told; fb = fold; dfb = dfold;
            exitflag = 2;
            return
        end

        % Update t
        if 2*t - told < tmax % if t + (t - told) < tmax
            lb = 2*t-told; % lb = t + (t - told) >= tmax
            ub = min(tmax,t+tau1*(t-told));
            tnew = FletcherCubInterp(told,fold,dfold,t,ft,dft);
            tnew = min(max(tnew,lb),ub);
            told = t;
            t = tnew;
        else
        	t = tmax;
        end
    end

    % We reach this point if and only if maxnf was reached
    exitflag = -1;
end

%-----------------------------------------------------------------------------------
function [cachet,t,exitflag,cnt] = FletcherSection(cache,prob,gam,f0,a,fa,dfa,b,fb,dfb,lsopt) 
%
% sectioningPhase finds an acceptable point t within a given bracket [a,b] 
% containing acceptable points. Notice that funcCount counts the total number of 
% function evaluations including those of the bracketing phase. 

    %      Q, A, C, f2, proxg
    cnt = [0, 0, 0, 0, 0];
    tau2 = min(0.1, lsopt.sigma); 
    tau3 = 0.5;

    t = [];
    for nsection = 1:lsopt.nsection 

        % Pick t in reduced interval
        lb = a + tau2*(b - a); ub = b - tau3*(b - a);
        % Find global minimizer in [lb, ub] of 3rd-degree polynomial that interpolates
        % f() and f'() at "a" and at "b".
        t = FletcherCubInterp(a,fa,dfa,b,fb,dfb);
        t = min(max(t,lb),ub);

        [cachet, cnt1] = DirFBE(prob, gam, t, cache, 3);
        cnt = cnt+cnt1;
        ft = cachet.FBE; dft = cachet.dFBE;

        if (t - a)*dfa >= lsopt.progTol || abs(b - a)*norm(cache.dir) < lsopt.progTol
          exitflag = -2;  % No further progress can be made
          return
        end

        % Update bracket
        aold = a; faold = fa; dfaold = dfa;  
        bold = b; fbold = fb; dfbold = dfb;
        if ft > f0 + t*lsopt.wolfe_hi || ft >= fa
            a = aold; fa = faold;  dfa = dfaold;
            b = t;  fb = ft; dfb = dft;
        else
            if abs(dft) <= -lsopt.wolfe_lo
                exitflag = 0; % Acceptable point found
                return
            end
            a = t; fa = ft; dfa = dft;
            if (b - a)*dft >= 0
                b = aold; fb = faold; dfb = dfaold;
            else
                b = bold; fb = fbold; dfb = dfbold;
            end
        end
    end % of while

    % We reach this point if and only if maxnf was reached
    exitflag = -1;
end

function t = FletcherCubInterp(t1,f1,df1,t2,f2,df2)
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

