%SEPARABLESUM Combines functions into their separable sum
%
%   SEPARABLESUM(fs, sizes, idx) where fs is a cell array of function
%   objects, sizes is a cell array containing size vectors, idx is an
%   integer vector of the same length of sizes.
%   If length(idx) = length(sizes) = k, then SEPARABLESUM returns the
%   function object correspondent to the sum
%
%       f(x) = sum_i=1...k fs{idx(i)}(x_i)
%
%   i.e., the sum of k functions, the ith being idx(i) and applied to a
%   block of prod(sizes{i}) variables.

% Copyright (C) 2015-2016, Lorenzo Stella and Panagiotis Patrinos
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

function obj = separableSum(objs, dims, idx)

    l = length(objs);
    if nargin < 3, idx = 1:l; end
    if length(idx) == 1, obj = objs{idx(1)}; return; end
    % determine Lipschitz constant (if possible)
    maxL = -1;
    noL = 0;
    obj.isConvex = 1;
    obj.isQuadratic = 1;
    obj.isConjQuadratic = 1;
    obj.hasHessian = 1;
    for i = 1:l
        if ~isfield(objs{i}, 'L'), noL = 1;
        elseif objs{i}.L > maxL, maxL = objs{i}.L; end
        if ~isfield(objs{i}, 'isConvex') || objs{i}.isConvex == 0
            obj.isConvex = 0;
        end
        if ~isfield(objs{i}, 'isQuadratic') || objs{i}.isQuadratic == 0
            obj.isQuadratic = 0;
        end
        if ~isfield(objs{i}, 'isConjQuadratic') || objs{i}.isConjQuadratic == 0
            obj.isConjQuadratic = 0;
        end
        if ~isfield(objs{i}, 'hasHessian') || objs{i}.hasHessian == 0
            obj.hasHessian = 0;
        end
    end
    if noL == 0
        obj.L = maxL;
    end
    for i = 1:length(dims)
        if numel(dims{i}) == 1, dims{i} = [dims{i}, 1]; end
    end
    nsum(1) = prod(dims{1});
    for i = 2:length(idx)
        nsum(i) = nsum(i-1) + prod(dims{i});
    end
    obj.makeprox = @() make_separableSum_prox(objs, idx, nsum, dims);
    obj.makef = @() make_separableSum_callf(objs, idx, nsum, dims);

end

function fun = make_separableSum_prox(objs, idx, nsum, dims)
    proxes = {};
    for i=1:length(objs)
        proxes{end+1} = objs{i}.makeprox();
    end
    fun = @(x, gam) call_separableSum_prox(x, gam, proxes, idx, nsum, dims);
end

function [prox, val] = call_separableSum_prox(x, gam, proxes, idx, nsum, dims)
    prox = zeros(nsum(end), 1);
    val = 0;
    baseidx = 1;
    for i=1:length(idx)
        xcurr = x(baseidx:nsum(i));
        [z, val1] = proxes{idx(i)}(reshape(xcurr, dims{i}(1), dims{i}(2)), gam);
        prox(baseidx:nsum(i)) = z(:);
        val = val+val1;
        baseidx = nsum(i)+1;
    end
end

function fun = make_separableSum_callf(objs, idx, nsum, dims)
    callfs = {};
    for i=1:length(objs)
        callfs{end+1} = objs{i}.makef();
    end
    fun = @(x) call_separableSum_f(x, callfs, idx, nsum, dims);
end

function [val, grad] = call_separableSum_f(x, callfs, idx, nsum, dims)
    val = 0;
    grad = [];
    baseidx = 1;
    if nargout >= 2
        for i=1:length(idx)
            xcurr = x(baseidx:nsum(i));
            [val1, grad1] = callfs{idx(i)}(reshape(xcurr, dims{i}(1), dims{i}(2)));
            val = val+val1;
            grad = [grad; grad1(:)];
            baseidx = nsum(i)+1;
        end
    else
        for i=1:length(idx)
            xcurr = x(baseidx:nsum(i));
            [val1] = callfs{idx(i)}(reshape(xcurr, dims{i}(1), dims{i}(2)));
            val = val+val1;
            baseidx = nsum(i)+1;
        end
    end
end
