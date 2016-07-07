% TODO: Problem structure, most general case.
%
% The structure of the problem is similar to the one discussed in the TFOCS
% user guide, available at http://cvxr.com/tfocs/doc/. (See in particular
% sec. 3 for the functions implemented there, and sec. 4.3 for the problem
% structure).
%
% Composite case:
%
%     f_1 (C_11 x_1 + ... + C1N x_N + d_1)
%   + f_2 (C_21 x_1 + ... + C2N x_N + d_2)
%   + ...
%   + f_M (C_M1 x_1 + ... + CMN x_N + d_M)
%   + g_1 (x_1) + ... + g_N (x_N)
%
% The terms in the cost must be split into quadratic/non-quadratic (the f_i
% or their conjugates) which will be then assigned to prob.f1 (quadratics)
% and to prob.f2 (non-quadratics). Accordingly, the C_ij and d_i will be
% stacked, and the g_i separable-summed with the appropriate dimensions.
%
% Separable case:
%
%   f_1 (x_1) + ... + f_M (x_M) + g_1 (z_1) + ... + g_N (z_N)
%
%   A_11 x_1 + ... + A_1M x_M + B_1 z_1 = b_1
%   A_21 x_1 + ... + A_2M x_M + B_2 z_2 = b_2
%   ...
%   A_N1 x_1 + ... + A_NM x_M + B_N z_N = b_N
%
% In this case the dual cost is computed as
%
%   min_{x_1} f_1(x_1) + y_1'(A_11 x_1) + ... + y_N'(A_N1 x_1)
%       = f*_1(-A*_11 y_1 - ... - A*_N1 y_N)
%
% (same for all f_i and g_i) so that the dual problem becomes
%
%     f*_1 (-A*_11 y_1 - ... - A*_N1 y_N)
%   + f*_2 (-A*_12 y_1 - ... - A*_N2 y_N)
%   + ...
%   + f*_M (-A*_1M y_1 - ... - A*_NM y_N)
%   + g*_1 (-B*_1 y_1) + ... + g*_N (-B*_N y_N)
%   + b_1'y_1 + ... + b_N'y_N
%
% which can be solved with forward-backward if B*_i B_i = multiple of the
% identity, for all i = 1,...,N.

function prob = MakeProblem(fs, gs, init, aff, constr)
    M = length(fs);
    N = length(gs);

    if M ~= 1
        error('only one smooth function is currently supported');
    end

    if N ~= 1
        error('only one nonsmooth function is currently supported');
    end

    if ~isa(fs, 'cell'), fs = {fs}; end
    if ~isa(gs, 'cell'), gs = {gs}; end

    if isempty(constr), flagconstr = 0; else flagconstr = 1; end
    if isempty(aff), flagaff = 0; flagd = 0; else flagaff = 1; end

    if flagconstr && flagaff
        error('cannot have both constraints and affine mappings');
    end

    if flagconstr
        if ~isa(constr, 'cell')
            error('the constraint must be a cell array');
        end
        if length(constr) ~= N+M+1
            error('must have as many terms in the constraints as f and g functions, plus the rhs');
        end
    end

    if flagaff
        if isa(aff, 'double') || isa(aff, 'struct')
            aff = {aff};
        end
        if ~isa(aff, 'cell')
            error('the list of affine maps must be a cell array or a matrix');
        end
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
    
    for i = 1:M, fs{i} = ProcessFunction(fs{i}); end
    for i = 1:N, gs{i} = ProcessFunction(gs{i}); end

    switch prob.id

        case 1 % composite problem (maybe aff, no constr)
            prob.x0 = init;
            for i = 1:M
                if fs{i}.isQuadratic
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
                else
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
            gs{1} = ProcessFunction(gs{1});
            prob.g = gs{1};

        case 2 % separable problem (constr, no aff)
            prob.y0 = init;
            for i = 1:M
                if fs{i}.isConjQuadratic
                    prob.f1 = fs{1};
                    if isstruct(constr{i})
                        prob.A1 = constr{i}.makeop();
                        prob.A1t = constr{i}.makeadj();
                    else
                        prob.A1 = constr{i};
                    end
                else
                    prob.f2 = fs{i};
                    if isstruct(constr{i})
                        prob.A2 = constr{i}.makeop();
                        prob.A2t = constr{i}.makeadj();
                    else
                        prob.A2 = constr{i};
                    end
                end
            end
            gs{1} = ProcessFunction(gs{1});
            prob.g = gs{1};
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