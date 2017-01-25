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
        % [obj.LRs, obj.Ks, obj.Ms, obj.Ls] = RiccatiFactor(obj.Q, obj.R, obj.Q_f, obj.A, obj.B, obj.N);
        obj = RiccatiFactor(obj);
        % obj.makefconj = @() make_lqrCost_fconj(x0, obj);
        obj.makefconj = @() @(w) call_lqrCost_fconj(w, x0, obj);
    else
        obj = varargin{1};
        % obj.makefconj = @() make_lqrCost_fconj(x0, obj);
        obj.makefconj = @() @(w) call_lqrCost_fconj(w, x0, obj);
    end
    obj.isConjQuadratic = 1;
end

function [fcw, xu] = call_lqrCost_fconj(w, x0, obj)
    [n_x, n_u] = size(obj.B);
    [~, xu] = RiccatiSolve(w, x0, obj.A, obj.B, obj.LRs, obj.Ks, obj.Ms, obj.Ls, int32(n_x), int32(n_u), int32(obj.N));
    fxu = 0;
    for i=0:obj.N-1
        x_i = xu(i*(n_x+n_u)+1:i*(n_x+n_u)+n_x);
        u_i = xu(i*(n_x+n_u)+n_x+1:(i+1)*(n_x+n_u));
        fxu = fxu + 0.5*(x_i'*(obj.Q*x_i) + u_i'*(obj.R*u_i));
    end
    x_N = xu(obj.N*(n_x+n_u)+1:end);
    fxu = fxu + 0.5*(x_N'*(obj.Q_f*x_N));
    fcw = w'*xu - fxu;
end

function obj = RiccatiFactor(obj)
    n = size(obj.Q,1);
    m = size(obj.R,1);
    Ps = zeros(n, n, obj.N+1);
    Ps(:,:,obj.N+1) = obj.Q_f;
    obj.LRs = zeros(m, m, obj.N);
    obj.Ss = zeros(m, n, obj.N);
    obj.Ks = zeros(m, n, obj.N);
    obj.Ms = zeros(m, n, obj.N);
    obj.Ls = zeros(n, n, obj.N);
    for k = obj.N:-1:1
        Rbar = obj.R + obj.B'*(Ps(:,:,k+1)*obj.B);
        Rbar = (Rbar+Rbar')/2;
        LR = chol(Rbar, 'lower');
        obj.LRs(:,:,k) = LR;
        obj.Ss(:,:,k) = obj.B'*(Ps(:,:,k+1)*obj.A);
        obj.Ks(:,:,k) = -(LR'\(LR\obj.Ss(:,:,k)));
        Ps(:,:,k) = obj.Q + obj.A'*(Ps(:,:,k+1)*obj.A) + obj.Ss(:,:,k)'*obj.Ks(:,:,k);
        Ps(:,:,k) = (Ps(:,:,k) + Ps(:,:,k)')/2;
    end
    for k = 1:obj.N
        LR = obj.LRs(:,:,k);
        obj.Ms(:,:,k) = -(LR'\(LR\obj.B'));
        obj.Ls(:,:,k) = (obj.A + obj.B*obj.Ks(:,:,k))';
    end
end
