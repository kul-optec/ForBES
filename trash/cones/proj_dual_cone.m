function x = proj_dual_cone(x, K)

% y = proj_cone(-x, K);
% x = x + y; % Moreau identity

% It looks like for some cones (such as second-order cones) project directly
% onto the dual cone when one knows how to do that, rather than using
% Moreau identity, yields more accurate results.

% dual zero cone (i.e. the free cone)
% nothing to do
idx = K.zero;

% nonnegative orthant (self-dual)
x(idx+1:idx+K.nn) = max(x(idx+1:idx+K.nn), 0);
idx = idx + K.nn;

% second-order cones (self-dual)
for i=1:length(K.soc)
    x(idx+1:idx+K.soc(i)) = proj_soc(x(idx+1:idx+K.soc(i)));
    idx = idx + K.soc(i);
end

end
