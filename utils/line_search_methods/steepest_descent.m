function [xmin, ymin] = steepest_descent(f, f_grad, x, line_search_method, a=1, rho=0.5, c=0.1, eps=1e-6, max_iters=100)
    # Line search using the steepest descent method

    fprintf("STARTED Line search using steepest descent\n")
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
        fprintf("Iter %d, x=[%f,%f], a=%f\n", iter, x(1), x(2), step)
    end
    xmin = x;
    ymin = f(xmin);
    fprintf("ENDED Line search using steepest descent\n")
end

