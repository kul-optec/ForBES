function [y, v] = prox(obj, x, gam)
    switch obj.mode
        case 'exact' % exact prox
            [U, S, V] = svd(x, 'econ');
            diagS1 = max(0, diag(S)-obj.lam*gam);
            S1 = diag(sparse(diagS1));
            y = U*(S1*V');
            v = obj.lam*sum(diagS1);
        otherwise
            % TODO: implement these
            error('not supported');
    end
end
