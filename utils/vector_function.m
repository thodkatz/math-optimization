function ret = vector_function(fun_handler, vector)
    if nargin(fun_handler) ~= 0
        ret = fun_handler(num2cell(vector){:});
    else
        ret = fun_handler();
    end
end