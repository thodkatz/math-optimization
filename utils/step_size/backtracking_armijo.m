function a = backtracking_armijo(f, f_grad, x, pk, a, rho, c, iter_count=10)
    % Backtracking line search using the sufficient decrease condition (armijo).

    iter = 0;
    % fprintf("\nSTARTED ARMIJO Step size config...\n")
    ck = f(x);
    phi0 = dot(f_grad(x), pk);
    while armijo_condition(f, phi0, x, pk, a, rho, c, ck)
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

function b = armijo_condition(f, phi0, x, pk, a, rho, c, ck)
    b = f(x + a*pk) > ck + c*a*phi0;
end