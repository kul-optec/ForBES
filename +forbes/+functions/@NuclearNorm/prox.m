function [y, v] = prox(obj, x, gam)
    switch obj.mode
        case 'exact' % exact prox
            [U, S, V] = svd(x, 'econ');
            diagS1 = max(0, diag(S)-obj.lam*gam);
            S1 = diag(sparse(diagS1));
            y = U*(S1*V');
            v = obj.lam*sum(diagS1);
        case 'inexact'
            [m, n] = size(x);
            maxrank = min(m, n);
            [U, S, V] = svds(x, obj.nsv, 'L');
            diagS1 = max(0, diag(S)-obj.lam*gam);
            if nnz(diagS1) == length(diagS1)
                obj.nsv = min(maxrank, obj.nsv+5);
            else
                obj.nsv = nnz(diagS1)+1;
            end
            S1 = diag(sparse(diagS1));
            y = U*(S1*V');
            if nargout >= 2
                v = obj.lam*sum(diagS1);
            end
        otherwise
            % TODO: implement these
            error('not supported');
    end
end
