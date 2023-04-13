function step = step_size(f, f_grad, line_searh_method, x, pk, a, rho, c)
    is_wolfe = 1;
    is_armijo = 2;
    is_none = 3;
    if line_searh_method == is_armijo
        step = backtracking_armijo(f, f_grad, x, pk, a, rho, c);
    elseif line_searh_method == is_wolfe
        step = wolfe(f, f_grad, x, pk, a, rho, c);
    else
        step = a;
    end
end