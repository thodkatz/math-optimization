function step = step_size(f, f_grad, line_search_method, x, pk, a, rho, c)
    is_backtracking_wolfe_weak = 1;
    is_wolfe_strong = 2;
    is_backtracking_armijo = 3;
    is_none = 4;
    if line_search_method == is_backtracking_armijo
        step = backtracking_armijo(f, f_grad, x, pk, a, rho, c);
    elseif line_search_method == is_backtracking_wolfe_weak
        step = backtracking_wolfe_weak(f, f_grad, x, pk, a, rho, c);
    elseif line_search_method == is_wolfe_strong
        step = wolfe_strong(f, f_grad, x, pk, a, rho, c);
    else
        step = a;
    end
end