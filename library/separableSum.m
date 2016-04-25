%SEPARABLESUM Combines separable functions into their sum
%
%   SEPARABLESUM(fs, sizes, idx) where fs is a cell array of function
%   objects, while idx and sizes are integer vectors of the same length.
%   If length(idx) = length(sizes) = k, then SEPARABLESUM returns the
%   function object correspondent to the sum
%       
%       f(x) = sum_i=1...k fs{idx(i)}(x_i)
%
%   i.e., the sum of k functions, the ith being idx(i) and applied to a
%   block of size(i) variables.
%
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

function obj = separableSum(objs, sizes, idx)
    l = length(objs);
    if nargin < 3
        idx = 1:l;
    end
    % determine Lipschitz constant (if possible)
    maxL = -1;
    noL = 0;
    obj.isConvex = 1;
    for i = 1:length(objs)
        if ~isfield(objs{i}, 'L'), noL = 1; break; end
        if objs{i}.L > maxL, maxL = objs{i}.L; end
            if ~isfield(objs{i}, 'isConvex') || objs{i}.isConvex == 0
                obj.isConvex = 0;
            end
    end
    if noL == 0
        obj.L = maxL;
    end
    obj.makeprox = @() make_separableSum_prox(objs, idx, sizes);
    obj.makef = @() make_separableSum_callf(objs, idx, sizes);
end

function op = make_separableSum_prox(objs, idx, sizes)
    proxes = {};
    for i=1:length(objs)
        proxes{end+1} = objs{i}.makeprox();
    end
    op = @(x, gam) call_separableSum_prox(x, gam, proxes, idx, sizes);
end

function [prox, val] = call_separableSum_prox(x, gam, proxes, idx, sizes)
    n = sum(sizes);
    prox = zeros(n, 1);
    val = 0;
    baseidx = 0;
    for i=1:length(idx)
        endcurr = baseidx+sizes(i);
        xcurr = x(baseidx+1:endcurr);
        [prox(baseidx+1:endcurr), val1] = proxes{idx(i)}(xcurr, gam);
        val = val+val1;
        baseidx = endcurr;
    end
end

function op = make_separableSum_callf(objs, idx, sizes)
    callfs = {};
    for i=1:length(objs)
        if isfield(objs{i}, 'isQuadratic') && objs{i}.isQuadratic
            callfs{end+1} = @(x) call_quadratic_f(x, objs{i}.Q, objs{i}.q);
        else
            callfs{end+1} = objs{i}.makef();
        end
    end
    op = @(x) call_separableSum_f(x, callfs, idx, sizes);
end

function [val, grad] = call_separableSum_f(x, callfs, idx, sizes)
    n = sum(sizes);
    prox = zeros(n, 1);
    val = 0;
    grad = [];
    baseidx = 0;
    if nargout > 1
        for i=1:length(idx)
            endcurr = baseidx+sizes(i);
            xcurr = x(baseidx+1:endcurr);
            [val1, grad1] = callfs{idx(i)}(xcurr);
            val = val+val1;
            grad = [grad; grad1];
            baseidx = endcurr;
        end
    else
        for i=1:length(idx)
            endcurr = baseidx+sizes(i);
            xcurr = x(baseidx+1:endcurr);
            val1 = callfs{idx(i)}(xcurr);
            val = val+val1;
            baseidx = endcurr;
        end
    end
end

function [val, grad] = call_quadratic_f(x, Q, q)
    grad = Q*x+q;
    val = 0.5*(x'*(grad+q));
end