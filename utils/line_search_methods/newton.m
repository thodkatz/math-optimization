function [xmin, ymin] = newton(f_sym, x, line_search_method, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100)
    # Line search using the newton method
    f = @(v) sym2fun(f_sym, v);
    f_grad = @(v) multidim_grad(f_sym, v);
    f_hessian = @(v) multidim_hessian(f_sym, v);

    fprintf("STARTED Line search using newton\n")
    for iter = 1:max_iters
        if any(isinf(x))
            error("Failed to converge. Inf value reached")
        end

        grad_x = f_grad(x);
        hess_x = f_hessian(x);
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