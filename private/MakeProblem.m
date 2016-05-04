function prob = MakeProblem(fs, gs, init, aff, constr, opt)
    M = length(fs);
    N = length(gs);
    if M > 2
        error('only the sum of two functions is currently supported');
    end
    if ~isa(fs, 'cell'), fs = {fs}; end
    if isempty(constr)
        flagconstr = 0;
    else
        if ~isa(constr, 'cell')
            error('the constraint must be a cell array');
        end
        flagconstr = 1;
        if length(constr) ~= N+M+1
            error('must have as many terms in the constraints as f and g functions, plus the rhs');
        end
    end
    if isempty(aff)
        flagaff = 0;
        flagd = 0;
        ns = length(init);
    else
        if isa(aff, 'double') || isa(aff, 'struct')
            aff = {aff};
        end
        if ~isa(aff, 'cell')
            error('the list of affine maps must be a cell array, a matrix or a structure');
        end
        flagaff = 1;
        if length(aff) == N, flagd = 0;
        elseif length(aff) == N+1, flagd = 1;
        else error('must have as many blocks of variables as g functions, if any');
        end
        ns = zeros(1,N);
        for i=1:N
            if ismatrix(aff{1,i}), ns(i) = size(aff{1,i},2);
            else ns(i) = aff{1,i}.n; end
        end
    end
    prob.id = flagconstr+1;
    switch prob.id
        case 1
            prob.x0 = init;
            for i = 1:M
                if isfield(fs{i}, 'isQuadratic') && fs{i}.isQuadratic
                    prob.f1 = fs{i};
                    if flagaff
                        op = separableSumLinear(aff(i,1:N));
                        if isstruct(op)
                            prob.C1 = op.makeop();
                            prob.C1t = op.makeadj();
                        else
                            prob.C1 = op;
                        end
                        if flagd, prob.d1 = -aff{i,N+1}; end
                    end
                end
                if ~isfield(fs{i}, 'isQuadratic') || ~fs{i}.isQuadratic
                    prob.f2 = fs{i};
                    if flagaff
                        op = separableSumLinear(aff(i,1:N));
                        if isstruct(op)
                            prob.C2 = op.makeop();
                            prob.C2t = op.makeadj();
                        else
                            prob.C2 = op;
                        end
                        if flagd, prob.d2 = -aff{i,N+1}; end
                    end
                end
            end
            if N == 1
                if isa(gs, 'cell'), prob.g = gs{1};
                else prob.g = gs; end
            elseif N > 1
                error('not implemented yet');
            end
        case 2
            prob.y0 = init;
            if flagaff
                error('cannot have both constraints and affine mappings');
            end
            for i = 1:M
                if isfield(fs{i}, 'isConjQuadratic') && fs{i}.isConjQuadratic
                    prob.f1 = fs{1};
                    if isstruct(constr{i})
                        prob.A1 = constr{i}.makeop();
                        prob.A1t = constr{i}.makeadj();
                    else
                        prob.A1 = constr{i};
                    end
                end
                if ~isfield(fs{i}, 'isConjQuadratic') || ~fs{i}.isConjQuadratic
                    prob.f2 = fs{i};
                    if isstruct(constr{i})
                        prob.A2 = constr{i}.makeop();
                        prob.A2t = constr{i}.makeadj();
                    else
                        prob.A2 = constr{i};
                    end
                end
            end
            if N == 1
                if isa(gs, 'cell'), prob.g = gs{1};
                else prob.g = gs; end
            elseif N > 1
                error('not implemented yet');
            end
            prob.B = horzcat(constr{M+1:M+N});
            prob.b = constr{M+N+1};
    end
end

function op = separableSumLinear(linops)
    N = length(linops);
    % for a single linear operator, return the operator itself
    if N == 1
        if isa(linops, 'cell'), op = linops{1};
        else op = linops; end
        return;
    end
    flagstruct = zeros(1,N);
    n = 0;
    for i = 1:N
        if isstruct(linops{i})
            flagstruct(i) = 1;
            m = linops{i}.m;
            n = n + linops{i}.n;
            linops{i}.op = linops{i}.makeop();
            linops{i}.adj = linops{i}.makeadj();
        else
            m = size(linops{i},1);
            n = n + size(linops{i},2);
        end
    end
    % if all operators are matrices, simply stack them horizontally
    if flagstruct == 0, op = horzcat(linops{:}); return; end
    % otherwise make a structure
    op.m = m;
    op.n = n;
    op.makeop = @() @(x) call_separableSumLinear(x, linops, N, m, flagstruct);
    op.makeadj = @() @(y) call_separableSumLinear_adj(y, linops, N, flagstruct);
end

function y = call_separableSumLinear(x, linops, N, m, flags)
    y = zeros(m,1);
    k = 0;
    for i=1:N
        if flags(i)
            y = y + linops{i}.op(x(k+1:k+linops{i}.n));
        else
            y = y + linops{i}*x(k+1:k+linops{i}.n);
        end
    end
end

function x = call_separableSumLinear_adj(y, linops, N, flags)
    x = [];
    for i=1:N
        if flags(i)
            x = [x; linops{i}.adj(y)];
        else
            x = [x; linops{i}*y];
        end
    end
end