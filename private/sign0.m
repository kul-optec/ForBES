% Custom sign function: sign0(x) = sign(x) if x != 0, 1 otherwise

function x = sign0(x)

if x == 0
    x = 1;
    return
end
x = sign(x);

end