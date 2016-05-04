%LQRCOST Allocates the linear quadratic regulator (LQR) cost function
%
%   LQRCOST(x0, Q, R, Q_f, A, B, N) builds the LQR cost with stage matrices
%   Q (for states) and R (for inputs), final cost matrix Q_f, dynamics A
%   and B, prediction horizon N and initial state x0.
%
%   LQRCOST(x0, obj) updates and return the LQR function obj with the new
%   initial state x0.
%
%   Example:
%       f = LQRCOST(x0, Q, R, Q_f, A, B, N);
%       [compute the next state x1 of the system]
%       f = LQRCOST(x1, f);
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

function obj = lqrCost(x0, varargin)
    %
    % Only f conjugate is available.
    %
    if length(varargin) > 1
        obj.Q = varargin{1};
        obj.R = varargin{2};
        obj.Q_f = varargin{3};
        obj.A = varargin{4};
        obj.B = varargin{5};
        obj.N = varargin{6};
        [obj.LRs, obj.Ks, obj.Ms, obj.Ls] = RiccatiFactor(obj.Q, obj.R, obj.Q_f, obj.A, obj.B, obj.N);
        obj.makefconj = @() make_lqrCost_fconj(x0, obj.A, obj.B, obj.N, obj.LRs, obj.Ks, obj.Ms, obj.Ls);
    else
        obj = varargin{1};
        obj.makefconj = @() make_lqrCost_fconj(x0, obj.A, obj.B, obj.N, obj.LRs, obj.Ks, obj.Ms, obj.Ls);
    end
    obj.isConjQuadratic = 1;
end

function op = make_lqrCost_fconj(x0, A, B, N, LRs, Ks, Ms, Ls)
    [n, m] = size(B);
    op = @(w) RiccatiSolve(w, x0, A, B, LRs, Ks, Ms, Ls, int32(n), int32(m), int32(N));
end

function [LRs, Ks, Ms, Ls] = RiccatiFactor(Q, R, Q_f, A, B, N)
    n = size(Q,1);
    m = size(R,1);
    Ps = zeros(n, n, N+1);
    Ps(:,:,N+1) = Q_f;
    LRs = zeros(m, m, N);
    Ss = zeros(m, n, N);
    Ks = zeros(m, n, N);
    Ms = zeros(m, n, N);
    Ls = zeros(n, n, N);
    for k = N:-1:1
        Rbar = R+B'*(Ps(:,:,k+1)*B);
        Rbar = (Rbar+Rbar')/2;
        LR = chol(Rbar, 'lower');
        LRs(:,:,k) = LR;
        Ss(:,:,k) = B'*(Ps(:,:,k+1)*A);
        Ks(:,:,k) = -(LR'\(LR\Ss(:,:,k)));
        Ps(:,:,k) = Q + A'*(Ps(:,:,k+1)*A) + Ss(:,:,k)'*Ks(:,:,k);
        Ps(:,:,k) = (Ps(:,:,k)+Ps(:,:,k)')/2;
    end
    for k = 1:N
        LR = LRs(:,:,k);
        Ms(:,:,k) = -(LR'\(LR\B'));
        Ls(:,:,k) = (A+B*Ks(:,:,k))';
    end
end
