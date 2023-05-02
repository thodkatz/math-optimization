function a = backtracking_wolfe_weak(f, f_grad, x, pk, a, rho, c=[1e-4, 0.9], iter_count=10)
    % Backtracking line search using the sufficient decrease condition and the curvature condition.

    c1 = c(1);
    c2 = c(2);

    iter = 0;
    % fprintf("\nSTARTED WOLFE Step size config...\n")
    ck = f(x);
    phi0 = dot(f_grad(x), pk);
    while armijo_condition(f, phi0, x, pk, a, rho, c1, ck) || curvature_condition(f_grad, x, pk, a, rho, c2, phi0)
        a *= rho;
        % fprintf("Iter %d, step size: %d\n", iter, a)
        iter += 1;
        if iter == iter_count
            % fprintf("Maximum number of iterations\n")
            break
        end
    end
    % fprintf("ENDED WOLFE Step size a=%f \n\n", a)
end

function b = curvature_condition(f_grad, x, pk, a, rho, c, phi0)
    b = dot(f_grad(x + a*pk), pk) < c * phi0;
end


function b = armijo_condition(f, phi0, x, pk, a, rho, c, ck)
    b = f(x + a*pk) > ck + c*a*phi0;
end