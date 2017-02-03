function x = proj_cone(x, K)

% zero cone
x(1:K.zero) = 0;
idx = K.zero;

% nonnegative orthant
x(idx+1:idx+K.nn) = max(x(idx+1:idx+K.nn), 0);
idx = idx + K.nn;

% second-order cones
for i=1:length(K.soc)
    x(idx+1:idx+K.soc(i)) = proj_soc(x(idx+1:idx+K.soc(i)));
    idx = idx + K.soc(i);
end

% % semidefinite cones
% for i=1:length(K.sdc)
%     z(idx+1:idx+getSdConeSize(K.sdc(i))) = proj_sdc(z(idx+1:idx+getSdConeSize(K.sdc(i))),K.sdc(i));
%     idx=idx+getSdConeSize(K.sdc(i));
% end

% % exponential cones
% for i = 1:K.exp
%     z(idx+1:idx+3) = proj_exp(z(idx+1:idx+3));
%     idx = idx+3;
% end

% % power cones
% for i = 1:length(K.pow)
%     if (K.pow(i) > 0)
%         % primal
%         z(idx+1:idx+3) = proj_pow(z(idx+1:idx+3), K.ppow(i));
%     else
%         % dual
%         z(idx+1:idx+3) = z(idx+1:idx+3) + proj_pow(-z(idx+1:idx+3), -K.pow(i));
%     end
%     idx = idx+3;
% end

end