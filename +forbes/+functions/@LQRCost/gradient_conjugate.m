function [xu, fcw] = gradient_conjugate(obj, w)
    [n_x, n_u] = size(obj.B);
    [~, xu] = forbes.utils.RiccatiSolve(w+obj.tilt, obj.x0, obj.A, obj.B, obj.LRs, obj.Ks, obj.Ms, obj.Ls, int32(n_x), int32(n_u), int32(obj.N));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Less efficient
%     fxu = 0;
%     for i=0:obj.N-1
%         x_i = xu(i*n_xu+1:i*n_xu+n_x);
%         u_i = xu(i*n_xu+n_x+1:(i+1)*n_xu);
%         fxu = fxu + 0.5*(x_i'*(obj.Q*x_i) + u_i'*(obj.R*u_i));
%     end
    % More efficient
    XU_stage = reshape(xu(1:end-n_x), n_x + n_u, obj.N);
    fxu = 0.5*sum(sum(XU_stage.*(obj.QR*XU_stage)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_N = xu(obj.N*(n_x+n_u)+1:end);
    fxu = fxu + 0.5*(x_N'*(obj.Q_f*x_N));
    fcw = (w+obj.tilt)'*xu - fxu - obj.diff;
end
