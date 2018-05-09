function [function_value,gradient] = rosen(initial_point)
    a=1;
    b=100;
    function_value =(a-initial_point(1))^2 + b*(initial_point(2)-initial_point(1))^2;

    if nargout > 1
       gradient = [-2*(a-(b+1)*initial_point(1)+b*initial_point(2)); 2*b*(initial_point(2)-initial_point(1)) ];
    end
end
