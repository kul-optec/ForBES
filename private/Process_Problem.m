function [prob, id] = Process_Problem(fs, gs, init, aff, constr)

    if ~isempty(aff) && ~isempty(constr)
        error('cannot have both constraints and affine mappings');
    end

    M = length(fs);
    N = length(gs);

    if ~isa(fs, 'cell'), fs = {fs}; end
    if ~isa(gs, 'cell'), gs = {gs}; end

    if ~isempty(aff)
        if isa(aff, 'double') || isa(aff, 'struct')
            aff = {aff};
        end
        if ~isa(aff, 'cell')
            error('the list of affine maps must be a cell array or a matrix');
        end
        if size(aff,2) ~= N+1
            error('affine term doesn''t match the number of g''s');
        end
        if size(aff,1) ~= M
            error('affine term doesn''t match the number of f''s');
        end
    end

    if isempty(constr)
        id = 1;
        [f1, C1, d1, f2, C2, d2, g] = combineTermsComposite(fs, gs, aff);
        prob = forbes.problems.ProblemComposite(f1, C1, d1, f2, C2, d2, g, [], [], init);
    else
        id = 2;
        if ~isa(constr, 'cell')
            error('the constraint must be a cell array');
        end
        if size(constr, 1) ~= N
            error('constraint doesn''t match the number of gs');
        end
        if size(constr, 2) ~= M+2
            error('constraint doesn''t match the number of fs');
        end
        [f1, A1, f2, A2, g, B, b] = combineTermsSeparable(fs, gs, constr);
        f1c = [];
        f2c = [];
        if ~isempty(f1), f1c = forbes.functions.Conjugate(f1); end
        if ~isempty(f2), f2c = forbes.functions.Conjugate(f2); end
        gc = forbes.functions.Conjugate(g);
        prob = forbes.problems.ProblemComposite(f1c, -A1', [], f2c, -A2', [], gc, -B', b, init);
    end

end

function [idx_quad, idx_nonquad] = splitSmooth(fs, conj)

    idx_quad = [];
    idx_nonquad = [];

    for i = 1:length(fs)
        if ~conj && fs{i}.is_quadratic()
            idx_quad(end+1) = i;
        elseif conj && fs{i}.is_generalized_quadratic() && fs{i}.is_strongly_convex()
            idx_quad(end+1) = i;
        else
            idx_nonquad(end+1) = i;
        end
    end

end

function [f1, C1, d1, f2, C2, d2, g] = combineTermsComposite(fs, gs, aff)

    f1 = []; C1 = []; d1 = [];
    f2 = []; C2 = []; d2 = [];
    g = [];

    M = length(fs);
    N = length(gs);

    if N > 1
        dims = {};
        for j=1:N, dims{j} = size(aff{1,j},2); end
        g = forbes.functions.SeparableSum(gs, dims);
    else
        g = gs{1};
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
            f1 = forbes.functions.SeparableSum(fs(idx_quad), dims);
        else
            f1 = fs{idx_quad(1)};
        end
        if ~isempty(aff1)
            C1 = vertcat(aff1{idx_quad,1});
            d1 = vertcat(aff1{idx_quad,2});
        end
    end
    if ~isempty(idx_nonquad)
        if length(idx_nonquad) > 1
            dims = {};
            for i=1:length(idx_nonquad), dims{i} = size(aff1{idx_nonquad(i),1},1); end
            f2 = forbes.functions.SeparableSum(fs(idx_nonquad), dims);
        else
            f2 = fs{idx_nonquad(1)};
        end
        if ~isempty(aff1)
            C2 = vertcat(aff1{idx_nonquad,1});
            d2 = vertcat(aff1{idx_nonquad,2});
        end
    end

end

function [f1, A1, f2, A2, g, B, b] = combineTermsSeparable(fs, gs, constr)

    f1 = []; A1 = [];
    f2 = []; A2 = [];
    g = []; B = [];
    b = [];

    M = length(fs);
    N = length(gs);

    if N > 1
        dims = {};
        for j=1:N, dims{j} = size(constr{j,M+1},2); end
        g = forbes.functions.SeparableSum(gs, dims);
    else
        g = gs{1};
    end
    B = blkdiag(constr{:,M+1});

    constr1 = {};
    for i=1:N
        constr1{i} = vertcat(constr{1:N,i});
    end

    [idx_quad, idx_nonquad] = splitSmooth(fs, 1);
    if ~isempty(idx_quad)
        if length(idx_quad) > 1
            dims = {};
            for i=1:length(idx_quad), dims{i} = size(constr1{idx_quad(i)},2); end
            f1 = forbes.functions.SeparableSum(fs(idx_quad), dims);
        else
            f1 = fs{idx_quad(1)};
        end
        A1 = horzcat(constr1{idx_quad});
    end
    if ~isempty(idx_nonquad)
        if length(idx_nonquad) > 1
            dims = {};
            for i=1:length(idx_nonquad), dims{i} = size(constr1{idx_nonquad(i)},2); end
            f2 = forbes.functions.SeparableSum(fs(idx_nonquad), dims);
        else
            f2 = fs{idx_nonquad(1)};
        end
        A2 = horzcat(constr1{idx_nonquad});
    end
    b = vertcat(constr{:,end});

end
