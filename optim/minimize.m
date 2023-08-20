function [xmin, fmin, iter, num_fun_evals, num_grad_fun_evals] = minimize(f_sym, domains, x, method, line_search_method, options, varargin)
%function [xmin, fmin, iter, num_fun_evals, num_grad_fun_evals] = newton(f_sym, domains, x, line_search_method, search_domain_x, search_domain_y, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=1000, to_plot=true, varargin)
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
    %         Possible values are: 'bisection_wolfe_weak', 2: 'wolfe_strong', 3: 'backtracking_armijo', and 4: 'none'
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
    NEWTON = 0;
    STEEPEST = 1;
    LBFGS = 2;

    if strcmp(method, "newton")
        method = NEWTON;
    elseif strcmp(method, "steepest")
        method = STEEPEST;
    elseif strcmp(method, "lbfgs")
        method = LBFGS;
    else
        error(-1)
    end


    f = function_handle(f_sym);
    f = @(v,n=0) vector_function(f,v, n);

    sym_vars = symvar(f_sym);

    f_grad = gradient(f_sym, sym_vars);
    f_grad = function_handle(f_grad);
    f_grad = @(v,n=0) vector_function(f_grad, v, n);

    if method == NEWTON
        f_hessian = hessian(f_sym, sym_vars);
        f_hessian = function_handle(f_hessian);
        f_hessian = @(v,n=0) vector_function(f_hessian, v, n);
    end

    num_fun_evals = 0; 
    num_grad_fun_evals = 0;

    if numel(varargin) == 4
        to_plot = true;
        assert(varargin{1},"search_x")
        search_domain_x = varargin{2};
        assert(varargin{3},"search_y")
        search_domain_y = varargin{4};
    elseif numel(varargin) ~= 0
        error("Unexpected number of variable arg inputs")
    elseif numel(varargin) == 0
        to_plot = false;
    end

    fprintf("STARTED Line search using newton\n")
    eps = options.eps;
    max_iters = options.max_iters;
    a = options.a;
    xall = {x};
    steps_all = [a];
    [fval, num_fun_evals] = f(x, num_fun_evals);
    fall = [fval];
    for iter = 1:max_iters
        if any(isinf(x))
            error("Failed to converge. Inf value reached")
        end

        [grad_x, num_grad_fun_evals] = f_grad(x, num_grad_fun_evals);
        if method == NEWTON
            [hess_x, num_grad_fun_evals] = f_hessian(x, num_grad_fun_evals);
            pk = -linsolve(hess_x, grad_x);
        elseif method == STEEPEST
            pk = -grad_x;
        else
            error(-1)
        end

        % ensure that the init alpha value will keep the x + alpha*pk, within the domain of the function and the derivative
        options.a = check_boundaries(x, a, pk, domains);
        

        % check for convergence
        if norm(grad_x) < eps
            break
        end

        if startsWith(line_search_method, "hanger_zhang")
            if iter == 1
                ck = fval;
                qk = 1;
            else
                [ck, qk] = hanger_zhang_attrs(fval, ck, qk); # full nonmonotone
            end
        elseif startsWith(line_search_method, "grippo");
            % source:
            % A NONMONOTONE LINE SEARCH TECHNIQUE
            % FOR NEWTON’S METHOD
            % L. GRIPPOf, F. LAMPARIELLOf AND S. LUCIDI"
            if iter == 1
                ck = fval;
            else
                % max of previous iterations
                ck = max(fall(max(1,end-options.memory_limit):end));
            end
        else
            ck = fval; # monotone
        end
        [step, num_fun_evals, num_grad_fun_evals] = step_size(f, f_grad, x, pk, line_search_method, options, ck, num_fun_evals, num_grad_fun_evals);
        x += step*pk;
        [fval, num_fun_evals] = f(x, num_fun_evals);
        fall(end+1) = fval;
        steps_all(end+1) = step;
        xall{end+1} = x;
        fprintf("Iter %d, x=[%0.4e,%0.4e], a=%0.4e, grad norm=%0.4e \n", iter, x(1), x(2), step, norm(pk))
    end
    fprintf("ENDED Line search using newton\n")
    xmin = x
    fmin = f(xmin)
    iter
    num_fun_evals
    num_grad_fun_evals

    if to_plot
        plot_line_search1(f, fall,xall, steps_all, search_domain_x, search_domain_y)
    end
end

