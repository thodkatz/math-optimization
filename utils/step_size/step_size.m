function [step, num_fun_evals, num_grad_fun_evals] = step_size(f, f_grad, x, pk, line_search_method, options, ck, num_fun_evals, num_grad_fun_evals)
    split_method = strsplit(line_search_method, "_");
    if startsWith(line_search_method, "hanger_zhang")
        line_search_method = strjoin(split_method(1,3:end), "_");
    elseif startsWith(line_search_method, "grippo")
        line_search_method = strjoin(split_method(1,2:end), "_");
    end

    if strcmp(line_search_method, "backtracking_armijo") == 1
        [step, num_fun_evals, num_grad_fun_evals] = backtracking_armijo(f, f_grad, x, pk, options, ck, num_fun_evals, num_grad_fun_evals);
    elseif strcmp(line_search_method, "bisection_wolfe_weak") == 1
        [step, num_fun_evals, num_grad_fun_evals] = bisection_curvature("weak wolfe", f, f_grad, x, pk, options, ck, num_fun_evals, num_grad_fun_evals);
    elseif strcmp(line_search_method, "bisection_goldstein") == 1
        [step, num_fun_evals, num_grad_fun_evals] = bisection_curvature("goldstein", f, f_grad, x, pk, options, ck, num_fun_evals, num_grad_fun_evals);
    elseif strcmp(line_search_method, "wolfe_strong") == 1
                [step, num_fun_evals, num_grad_fun_evals] = wolfe_strong(f, f_grad, x, pk, options, ck, num_fun_evals, num_grad_fun_evals);
    elseif strcmp(line_search_method, "none") == 1
        step = options.a;
    else
        error(-1)
    end
end
