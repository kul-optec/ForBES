function z = proj_soc(z)

if isempty(z)
    z=[];
    return;
elseif length(z)==1
    z = max(z,0);
    return;
end

v1 = z(1);
v2 = z(2:end);
normv2 = norm(v2);

if v1 <= -normv2
    z = zeros(length(z), 1);
elseif v1 >= normv2
    z = z;
else
    a = (v1+normv2)/2;
    z(1) = a;
    z(2:end) = a*(z(2:end)/normv2);
end

end