function a = backtracking_armijo(f, f_grad, x, pk, a, rho, c, iter_count=10)
    % Backtracking line search using the sufficient decrease condition (armijo).

    iter = 0;
    % fprintf("\nSTARTED ARMIJO Step size config...\n")
    while armijo_condition(f, f_grad, x, pk, a, rho, c)
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

function b = armijo_condition(f, f_grad, x, pk, a, rho, c)
    b = f(x + a*pk) > f(x) + c*a*dot(f_grad(x), pk);
end