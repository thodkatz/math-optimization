function a = backtracking_armijo(f, f_grad, x, pk, a, rho, c, ck, iter_count=10)
    % Jorge Nocedal and Stephen J Wright. Numerical optimization. Springer, 1999
    % page: 37
    %
    % Backtracking line search using the sufficient decrease condition (armijo).

    iter = 0;
    % fprintf("\nSTARTED ARMIJO Step size config...\n")
    phi_derivative0 = dot(f_grad(x), pk);
    while armijo_condition(f, phi_derivative0, x, pk, a, c, ck)
        a *= rho;
        % fprintf("Iter %d, step size: %d\n", iter, a)
        iter += 1;
        if iter == iter_count
            % fprintf("Maximum number of iterations\n")
            break
        end
    end
    % fprintf("ENDED ARMIJO Step size a=%f \n\n", a)
end

