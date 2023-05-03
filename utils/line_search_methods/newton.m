function [xmin, fmin] = newton(f_sym, x, line_search_method, search_domain_x, search_domain_y, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100)
    % NEWTON Find the minimum of a function using the Newton method
    %
    % USAGE:
    %       [xmin, fmin] = newton(f_sym, x, line_search_method, c, rho, a, eps, max_iters)
    %       finds a local minimum of the symbolic function f_sym, starting at point x,
    %       using the Newton method with the specified line search method.
    %
    % INPUTS:
    %       - f_sym: symbolic function of n variables, where n is the number of
    %         elements in x
    %       - x: initial point of the search (1-by-n row vector)
    %       - line_search_method: string indicating the line search method to use.
    %         Possible values are: 'backtracking_wolfe_weak', 2: 'wolfe_strong', 3: 'backtracking_armijo', and 4: 'none'
    %       - c: parameter for the Armijo and Wolfe search methods
    %       - rho: parameter for the Armijo and Wolfe line search methods
    %       - a: initial step size for the line search methods
    %       - eps: tolerance for the stopping criterion
    %       - max_iters: maximum number of iterations
    %
    % OUTPUTS:
    %       - xmin: 1-by-n row vector representing the point where the minimum is found
    %       - fmin: scalar value representing the minimum value found
    %
    % DESCRIPTION:
    %       The Newton method is a second-order optimization method that uses
    %       the Hessian matrix to guide the search towards the minimum. At each
    %       iteration, the method computes the Newton direction by solving the
    %       system of linear equations Hx = -g, where H is the Hessian matrix of
    %       f at the current point, and g is the gradient of f at that point.
    %       Then, a step size is computed using a line search method that
    %       satisfies the Wolfe conditions.
    %
    % EXAMPLES:
    %       syms x y
    %       f_sym = x^2 + y^2 + x*y + x + y + 1;
    %       [xmin, fmin] = newton(f_sym, [-1, -1], "backtracking_armijo", 'rho', 0.8, 'c', 0.2)
    %       % Output: xmin = [-0.4992, -0.4992], fmin = 0.7508
    %

    f = @(v) sym2fun(f_sym, v);
    f_grad = @(v) multidim_grad(f_sym, v);
    f_hessian = @(v) multidim_hessian(f_sym, v);

    fprintf("STARTED Line search using newton\n")
    xall = {x};
    steps_all = [a];
    for iter = 1:max_iters
        if any(isinf(x))
            error("Failed to converge. Inf value reached")
        end

        grad_x = f_grad(x);
        hess_x = f_hessian(x);
        pk = -linsolve(hess_x, grad_x);

        % check for convergence
        if norm(grad_x) < eps
            break
        end

        if startsWith(line_search_method, "nonmonotone")
            if iter == 1
                ck = f(x);
                qk = 1;
            else
                [ck, qk] = nonmonotone_attrs(f(x), ck, qk);
            end
        else
            ck = f(x);
        end
        step = step_size(f, f_grad, line_search_method, x, pk, a, rho, c, ck);
        x += step*pk;
        steps_all(end+1) = step;
        xall{end+1} = x;
        % fprintf("Iter %d, x=[%f,%f], a=%f\n", iter, x(1), x(2), step)
    end
    fprintf("ENDED Line search using newton\n")
    xmin = x
    fmin = f(xmin)

    plot_line_search1(f, xall, steps_all, search_domain_x, search_domain_y)
end