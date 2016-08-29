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

function prob = MakeProblem(fs, gs, init, aff, constr)

    if ~isempty(aff) && ~isempty(constr)
        error('cannot have both constraints and affine mappings');
    end

    M = length(fs);
    N = length(gs);

    if ~isa(fs, 'cell'), fs = {fs}; end
    if ~isa(gs, 'cell'), gs = {gs}; end

    for i = 1:M, fs{i} = ProcessFunction(fs{i}); end
    for i = 1:N, gs{i} = ProcessFunction(gs{i}); end

    if ~isempty(aff)
        if isa(aff, 'double') || isa(aff, 'struct')
            aff = {aff};
        end
        if ~isa(aff, 'cell')
            error('the list of affine maps must be a cell array or a matrix');
        end
        if size(aff,2) ~= N+1
            error('affine term doesn''t match the number of gs');
        end
        if size(aff,1) ~= M
            error('affine term doesn''t match the number of fs');
        end
    end

    if isempty(constr)
        prob = combineTermsComposite(fs, gs, aff);
        prob.x0 = init;
        prob.id = 1;
    else
        if ~isa(constr, 'cell')
            error('the constraint must be a cell array');
        end
        if size(constr, 1) ~= N
            error('constraint doesn''t match the number of gs');
        end
        if size(constr, 2) ~= M+2
            error('constraint doesn''t match the number of fs');
        end
        prob = combineTermsSeparable(fs, gs, constr);
        prob.y0 = init;
        prob.id = 2;
    end

end

function [idx_quad, idx_nonquad] = splitSmooth(fs, conj)

    idx_quad = [];
    idx_nonquad = [];

    for i = 1:length(fs)
        if ~conj && isfield(fs{i}, 'isQuadratic') && fs{i}.isQuadratic
            idx_quad(end+1) = i;
        elseif conj && isfield(fs{i}, 'isConjQuadratic') && fs{i}.isConjQuadratic
            idx_quad(end+1) = i;
        else
            idx_nonquad(end+1) = i;
        end
    end

end

function prob = combineTermsComposite(fs, gs, aff)

    M = length(fs);
    N = length(gs);

    if N > 1
        dims = {};
        for j=1:N, dims{j} = size(aff{1,j},2); end
        prob.g = separableSum(gs, dims);
    else
        prob.g = gs{1};
    end

    aff1 = {};
    if ~isempty(aff)
        for i=1:M
            aff1{i,1} = horzcat(aff{i,1:N});
            if length(aff(i,:)) == N, aff1{i,2} = 0;
            else aff1{i,2} = aff{i,N+1}; end
        end
    end

    [idx_quad, idx_nonquad] = splitSmooth(fs, 0);
    if ~isempty(idx_quad)
        if length(idx_quad) > 1
            dims = {};
            for i=1:length(idx_quad), dims{i} = size(aff1{idx_quad(i),1},1); end
            prob.f1 = separableSum(fs(idx_quad), dims);
        else
            prob.f1 = fs{idx_quad(1)};
        end
        if ~isempty(aff1)
            prob.C1 = vertcat(aff1{idx_quad,1});
            prob.d1 = vertcat(aff1{idx_quad,2});
        end
    end
    if ~isempty(idx_nonquad)
        if length(idx_nonquad) > 1
            dims = {};
            for i=1:length(idx_nonquad), dims{i} = size(aff1{idx_nonquad(i),1},1); end
            prob.f2 = separableSum(fs(idx_nonquad), dims);
        else
            prob.f2 = fs{idx_nonquad(1)};
        end
        if ~isempty(aff1)
            prob.C2 = vertcat(aff1{idx_nonquad,1});
            prob.d2 = vertcat(aff1{idx_nonquad,2});
        end
    end

end

function prob = combineTermsSeparable(fs, gs, constr)

    M = length(fs);
    N = length(gs);

    dims = {};
    for j=1:N, dims{j} = size(constr{j,M+1},2); end
    prob.g = separableSum(gs, dims);
    prob.B = blkdiag(constr{:,M+1});

    constr1 = {};
    for i=1:N
        constr1{i} = vertcat(constr{1:N,i});
    end

    [idx_quad, idx_nonquad] = splitSmooth(fs, 1);
    if ~isempty(idx_quad)
        dims = {};
        for i=1:length(idx_quad), dims{i} = size(constr1{idx_quad(i)},2); end
        prob.f1 = separableSum(fs(idx_quad), dims);
        prob.A1 = horzcat(constr1{idx_quad});
    end
    if ~isempty(idx_nonquad)
        dims = {};
        for i=1:length(idx_nonquad), dims{i} = size(constr1{idx_nonquad(i)},2); end
        prob.f2 = separableSum(fs(idx_nonquad), dims);
        prob.A2 = horzcat(constr1{idx_nonquad});
    end
    prob.b = vertcat(constr{:,end});

end
