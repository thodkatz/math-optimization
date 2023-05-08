function a = bisection_wolfe_weak(f, f_grad, x, pk, a, c=[1e-4, 0.9], ck, iter_count=10)
    % Implementation of bisection wolfe weak condition
    %
    % Convergence of descent methods with backtracking 
    % (Armijo) linesearch. Bisection algorithm for weak Wolfe
    % conditions, Anton Evgrafov

    c1 = c(1);
    c2 = c(2);

    amin = 0;
    amax = inf;

    iter = 0;
    % fprintf("\nSTARTED WOLFE Step size config...\n")
    phi_derivative0 = dot(f_grad(x), pk);
    while 1
        if armijo_condition(f, phi_derivative0, x, pk, a, c1, ck) 
            amax = a;
            a = (amin + amax)/2;
        elseif curvature_condition(f_grad, x, pk, a, c2, phi_derivative0)
            amin = a;
            if amax == Inf;
                a *= 2;
            else
                a = (amin+amax)/2;
            end
        else
            break
        end
    end
    % fprintf("ENDED WOLFE Step size a=%f \n\n", a)
end

function b = curvature_condition(f_grad, x, pk, a, c, phi0)
    b = dot(f_grad(x + a*pk), pk) < c * phi0;
end
