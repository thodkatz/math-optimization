function a = wolfe(f, f_grad, x, pk, a, rho, c1, c2=0.2, iter_count=20)
    % Backtracking line search using the sufficient decrease condition and the curvature condition.

    iter = 0;
    fprintf("\nSTARTED WOLFE Step size config...\n")
    while armijo_condition(f, f_grad, x, pk, a, rho, c1) || curvature_condition(f_grad, x, pk, a, rho, c2)
        a *= rho;
        fprintf("Iter %d, step size: %d\n", iter, a)
        iter += 1;
        if iter == iter_count
            fprintf("Maximum number of iterations\n")
            break
        end
    end
    fprintf("ENDED WOLFE Step size config...\n\n")
end

function b = curvature_condition(f_grad, x, pk, a, rho, c)
    b = dot(f_grad(x + a*pk), pk) < c * dot(f_grad(x), pk);
end


function b = armijo_condition(f, f_grad, x, pk, a, rho, c)
    b = f(x + a*pk) > f(x) + c*a*dot(f_grad(x), pk);
end