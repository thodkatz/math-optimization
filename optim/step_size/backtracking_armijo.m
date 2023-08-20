function [a, num_fun_evals, num_grad_fun_evals] = backtracking_armijo(f, f_grad, x, pk, options, ck, num_fun_evals, num_grad_fun_evals)
    % Jorge Nocedal and Stephen J Wright. Numerical optimization. Springer, 1999
    % page: 37
    %
    % Backtracking line search using the sufficient decrease condition (armijo).
    a = options.a;
    rho = options.rho;
    assert(rho<1)
    c = options.c;
    iter_count = options.max_iters_step_size;

    iter = 0;
    fprintf("\nSTARTED ARMIJO Step size config...\n")
    [fgrad, num_grad_fun_evals] = f_grad(x, num_grad_fun_evals);
    phi_derivative0 = dot(fgrad, pk);
    [condition, num_fun_evals] = armijo_condition(f, phi_derivative0, x, pk, a, c, ck, num_fun_evals);
    while condition
        a *= rho;
        fprintf("Iter %d, step size: %0.4e\n", iter, a)
        iter += 1;
        if iter == iter_count
            fprintf("Maximum number of iterations\n")
            break
        end
        [condition, num_fun_evals] = armijo_condition(f, phi_derivative0, x, pk, a, c, ck, num_fun_evals);
    end
    fprintf("ENDED ARMIJO Step size a=%f \n\n", a)
end

