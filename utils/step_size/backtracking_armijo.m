function [a, number_of_function_evaluations, number_of_gradient_function_evaluations] = backtracking_armijo(f, f_grad, x, pk, a, rho, c, ck, iter_count=10, number_of_function_evaluations, number_of_gradient_function_evaluations)
    % Jorge Nocedal and Stephen J Wright. Numerical optimization. Springer, 1999
    % page: 37
    %
    % Backtracking line search using the sufficient decrease condition (armijo).

    iter = 0;
    % fprintf("\nSTARTED ARMIJO Step size config...\n")
    [fgrad, number_of_gradient_function_evaluations] = f_grad(x, number_of_gradient_function_evaluations);
    phi_derivative0 = dot(fgrad, pk);
    [condition, number_of_function_evaluations] = armijo_condition(f, phi_derivative0, x, pk, a, c, ck, number_of_function_evaluations);
    while condition
        a *= rho;
        % fprintf("Iter %d, step size: %d\n", iter, a)
        iter += 1;
        if iter == iter_count
            % fprintf("Maximum number of iterations\n")
            break
        end
        [condition, number_of_function_evaluations] = armijo_condition(f, phi_derivative0, x, pk, a, c, ck, number_of_function_evaluations);
    end
    % fprintf("ENDED ARMIJO Step size a=%f \n\n", a)
end

