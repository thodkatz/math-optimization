function [xmin, ymin] = steepest_descent(f_sym, x, line_search_method, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100)
    # Line search using the steepest descent method
    f = @(v) sym2fun(f_sym, v);
    f_grad = @(v) multidim_grad(f_sym, v);

    fprintf("STARTED Line search using steepest descent\n")
    xall = {x};
    steps_all = [a];
    for iter = 1:max_iters
        if any(isinf(x))
            error("Failed to converge. Inf value reached")
        end

        pk = -f_grad(x);
        
        # check for convergence
        if norm(pk) < eps
            break
        end

        step = step_size(f, f_grad, line_search_method, x, pk, a, rho, c);
        x += step*pk;
        steps_all(end+1) = step;
        xall{end+1} = x;
        % fprintf("Iter %d, x=[%f,%f], a=%f\n", iter, x(1), x(2), step)
    end
    xmin = x;
    ymin = f(xmin);
    fprintf("ENDED Line search using steepest descent\n")

    plot_line_search(f, xall, steps_all)
end

