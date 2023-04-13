function [xmin, ymin] = newton(f, f_grad, x, line_search_method, a=1, rho=0.5, c=0.1, eps=1e-6, max_iters=10)
    # Line search using the newton method

    fprintf("STARTED Line search using newton\n")
    for iter = 1:max_iters
        if any(isinf(x))
            error("Failed to converge. Inf value reached")
        end

        grad_x = rosen_grad(x);
        hess_x = rosen_hessian(x);
        pk = -linsolve(hess_x, grad_x);

        # check for convergence
        if norm(grad_x) < eps
            break
        end

        step = step_size(f, f_grad, line_search_method, x, pk, a, rho, c);
        x += step*pk;
        fprintf("Iter %d, x=[%f,%f], a=%f\n", iter, x(1), x(2), step)
    end
    xmin = x;
    ymin = f(xmin);
    fprintf("ENDED Line search using newton\n")
end