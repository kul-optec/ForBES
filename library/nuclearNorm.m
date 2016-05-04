%NUCLEARNORM Allocates the nuclear norm function
%
%   NUCLEARNORM(m, n, lam, mode) builds the function
%       
%       g(x) = lam*||x||_*
%
%   where ||.||_* is the nuclear norm for m-by-n matrices, and x is assumed
%   to be a vector of length m*n, containing the stacked columns of an
%   m-by-n matrix. If the third argument lam is not provided, lam = 1.
%
%   Fourth argument 'mode' selects how to compute the proximal operator
%   associated with the function:
%    - mode == 'exact': compute the full svd using MATLAB's svd
%    - mode == 'adaptive': compute the exact prox using a partial svd
%    - mode == 'inexact': compute an inexact prox using a partial svd
%
%   Fifth argument 'method' allows to choose between MATLAB's svds and
%   PROPACK's lansvd when mode is 'adaptive' or 'inexact'
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

function obj = nuclearNorm(m, n, lam, mode, method)
    global nsv;
    global flagadd;
    nsv = 10;
    flagadd = 0;
    if nargin < 2
        error('you must provide the number of rows and columns, m and n, as arguments');
    end
    if nargin < 3, lam = 1; end
    if nargin < 4, mode = 'exact'; end
    if nargin < 5, method = 'svds'; end
    switch mode
        case 'exact' % exact prox
            obj.makeprox = @() @(x, gam) call_nuclearNorm_prox(x, gam, m, n, lam);
        case 'adaptive' % adaptive prox
            switch method
                case 'svds'
                    obj.makeprox = @() @(x, gam) call_nuclearNorm_prox_adaptive_svds(x, gam, m, n, lam);
                case 'lansvd'
                    obj.makeprox = @() @(x, gam) call_nuclearNorm_prox_adaptive_lansvd(x, gam, m, n, lam);
                otherwise
                    error('unknown method for computing SVDs');
            end
        case 'inexact' % adaptive/inexact prox
            switch method
                case 'svds'
                    obj.makeprox = @() @(x, gam) call_nuclearNorm_prox_inexact_svds(x, gam, m, n, lam);
                case 'lansvd'
                    obj.makeprox = @() @(x, gam) call_nuclearNorm_prox_inexact_lansvd(x, gam, m, n, lam);
                otherwise
                    error('unknown method for computing SVDs');
            end
    end
end

function [prox, val] = call_nuclearNorm_prox(x, gam, m, n, lam)
    [U, S, V] = svd(reshape(x, m, n), 'econ');
    diagS1 = max(0, diag(S)-lam*gam);
    S1 = diag(sparse(diagS1));
    prox = reshape(U*(S1*V'), m*n, 1);
    if nargout >= 2
        val = lam*sum(diagS1);
    end
end

function [prox, val] = call_nuclearNorm_prox_adaptive_svds(x, gam, m, n, lam)
    global nsv;
    global flagadd;
    maxrank = min(m, n);
    flagok = 0;
    while ~flagok
        [U, S, V] = svds(reshape(x, m, n), nsv, 'L');
        diagS1 = max(0, diag(S)-lam*gam);
        if nnz(diagS1) == length(diagS1)
            if flagadd
                nsv = min(maxrank, 10+nsv);
            else
                nsv = min(maxrank, 2*nsv);
            end
        else
            nsv = nnz(diagS1)+1;
            flagok = 1;
            flagadd = 1;
        end
    end
    S1 = diag(sparse(diagS1));
    prox = reshape(U*(S1*V'), m*n, 1);
    if nargout >= 2
        val = lam*sum(diagS1);
    end
end

function [prox, val] = call_nuclearNorm_prox_inexact_svds(x, gam, m, n, lam)
    global nsv;
    maxrank = min(m, n);
    [U, S, V] = svds(reshape(x, m, n), nsv, 'L');
    diagS1 = max(0, diag(S)-lam*gam);
    if nnz(diagS1) == length(diagS1)
        nsv = min(maxrank, nsv+5);
    else
        nsv = nnz(diagS1)+1;
    end
    S1 = diag(sparse(diagS1));
    prox = reshape(U*(S1*V'), m*n, 1);
    if nargout >= 2
        val = lam*sum(diagS1);
    end
end

function [prox, val] = call_nuclearNorm_prox_adaptive_lansvd(x, gam, m, n, lam)
    global nsv;
    global flagadd;
    maxrank = min(m, n);
    flagok = 0;
    while ~flagok
        [U, S, V] = lansvd(reshape(x, m, n), nsv, 'L');
        diagS1 = max(0, diag(S)-lam*gam);
        if nnz(diagS1) == length(diagS1)
            if flagadd
                nsv = min(maxrank, 10+nsv);
            else
                nsv = min(maxrank, 2*nsv);
            end
        else
            nsv = nnz(diagS1)+1;
            flagok = 1;
            flagadd = 1;
        end
    end
    S1 = diag(sparse(diagS1));
    prox = reshape(U*(S1*V'), m*n, 1);
    if nargout >= 2
        val = lam*sum(diagS1);
    end
end

function [prox, val] = call_nuclearNorm_prox_inexact_lansvd(x, gam, m, n, lam)
    global nsv;
    maxrank = min(m, n);
    [U, S, V] = lansvd(reshape(x, m, n), nsv, 'L');
    diagS1 = max(0, diag(S)-lam*gam);
    if nnz(diagS1) == length(diagS1)
        nsv = min(maxrank, nsv+5);
    else
        nsv = nnz(diagS1)+1;
    end
    S1 = diag(sparse(diagS1));
    prox = reshape(U*(S1*V'), m*n, 1);
    if nargout >= 2
        val = lam*sum(diagS1);
    end
end
