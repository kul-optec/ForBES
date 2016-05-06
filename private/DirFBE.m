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

function [cachet, ops] = DirFBE(cache, tau, mode, cachet)
% Computes the 'directional' FBE, i.e., FBE(x+tau*d) and its derivative with
% respect to tau, if requested. Here x = cache.x, d = cache.dir.
% If cachet (4th argument) is provided, then skips precomputing the data
% that has already been stored in cachet.
%
% If mode == 1, then compute only FBE(x+tau*d) and put it into cachet.
% If mode == 2, compute only dFBE(x+tau*d), the directional derivative.
% If mode == 3, compute both FBE and dFBE at x+tau*d.

if cache.flagLineSearch ~= 1
    error('DirFBE: line search data was not precomputed');
end

ops = OpsInit();

prob = cache.prob;
gam = cache.gam;

if nargin < 4
    cachet = CacheInit(prob, cache.x + tau*cache.dir, gam);
end

if nargin < 4 || cachet.flagGradStep ~= 1
    fxt = 0;
    gradfxt = 0;
    if prob.istheref1
        cachet.res1x = cache.res1x + tau*cache.C1dir;
        cachet.Qres1x = cache.Qres1x + tau*cache.QC1dir;
        cachet.gradf1x = cache.gradf1x + tau*cache.C1tQC1dir;
        cachet.f1x = cache.f1x + tau*cache.f1linear + (0.5*tau^2)*cache.f1quad;
        fxt = fxt + cachet.f1x;
        gradfxt = gradfxt + cachet.gradf1x;
    end
    if prob.istheref2
        cachet.res2x = cache.res2x + tau*cache.C2dir;
        if prob.useHessian
            [f2xt, gradf2res2xt, cachet.Hessf2res2x] = prob.callf2(cachet.res2x);
        else
            [f2xt, gradf2res2xt] = prob.callf2(cachet.res2x);
            cachet.gradf2res2x = gradf2res2xt;
        end
        ops.f2 = ops.f2 + 1;
        ops.gradf2 = ops.gradf2 + 1;
        if prob.isthereC2
            if prob.isC2fun, gradf2xt = prob.C2t(gradf2res2xt);
            else gradf2xt = prob.C2'*gradf2res2xt; end
            ops.C2 = ops.C2 + 1;
        else
            gradf2xt = gradf2res2xt;
        end
        fxt = fxt + f2xt;
        gradfxt = gradfxt + gradf2xt;
    end
    if prob.istherelin
        cachet.flinx = cache.flinx + tau*cache.lindir;
        fxt = fxt + cachet.flinx;
        gradfxt = gradfxt + prob.l;
    end
    % compute proximal gradient step
    cachet.fx = fxt;
    cachet.gradfx = gradfxt;
    cachet.y = cachet.x - gam*gradfxt;
    
    cachet.flagGradStep = 1;
end

if nargin < 4 || cachet.flagProxGradStep ~= 1
    [cachet.z, cachet.gz] = prob.callg(cachet.y, gam);
    ops.proxg = ops.proxg + 1;
    ops.g = ops.g + 1;
    cachet.diff = cachet.z-cachet.x;
    
    cachet.flagProxGradStep = 1;
end

if mode == 1 || mode == 3
    sqnormdifft = cachet.diff'*cachet.diff;
    cachet.normdiff = sqrt(sqnormdifft);
    cachet.FBE = cachet.fx + cachet.gz + cachet.gradfx'*cachet.diff + (0.5/gam)*sqnormdifft;
    
    cachet.flagFBE = 1;
end

if mode >= 2
    Hdir = 0;
    if prob.istheref1
        Hdir = Hdir + cache.C1tQC1dir;
    end
    if prob.istheref2
        if prob.useHessian
            HC2dir = cachet.Hessf2res2x*cache.C2dir;
        else
            res2xtepsdir = cachet.res2x + 1e-100i*cache.C2dir;
            [~, gradf2res2xtepsdir] = prob.callf2(res2xtepsdir);
            ops.f2 = ops.f2 + 1;
            HC2dir = imag(gradf2res2xtepsdir)/1e-100;
%                 res2xtepsdir = cache.res2x + 1e-8*cache.C2dir;
%                 [~, gradf2res2xtepsdir] = prob.callf2(res2xtepsdir);
%                 cnt(4) = cnt(4)+1;
%                 HC2dir = (gradf2res2xtepsdir-cachet.gradf2res2x)/1e-8;
        end
        if prob.isthereC2
            if prob.isC2fun, Hdir = Hdir + prob.C2t(HC2dir);
            else Hdir = Hdir + (prob.C2'*HC2dir); end
            ops.C2 = ops.C2 + 1;
        else
            Hdir = Hdir + HC2dir;
        end
    end
    cachet.dFBE = cachet.diff'*Hdir-(cachet.diff'*cache.dir)/gam;
end
