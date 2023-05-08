function step = step_size(f, f_grad, line_search_method, x, pk, a, rho, c, ck)
    split_method = strsplit(line_search_method, "_");
    if strcmp(split_method(1,1), "nonmonotone")
        line_search_method = strjoin(split_method(1,2:end), "_");
    end

    if strcmp(line_search_method, "backtracking_armijo") == 1
        step = backtracking_armijo(f, f_grad, x, pk, a, rho, c, ck);
    elseif strcmp(line_search_method, "bisection_wolfe_weak") == 1
        if !is_vector(c)
            c = [1e-4, 0.9];
        end
        step = bisection_wolfe_weak(f, f_grad, x, pk, a, c, ck);
    elseif strcmp(line_search_method, "wolfe_strong") == 1
        if !is_vector(c)
            c = [1e-4, 0.9];
        end
        if rho < 1
            fprintf("WARNING: rho < 1, should be greater for a to reach the amax\n")
            rho = 2;
        end
        step = wolfe_strong(f, f_grad, x, pk, a, rho, c, ck);
    elseif strcmp(line_search_method, "none") == 1
        step = a;
    else
        error(-1)
    end
end
