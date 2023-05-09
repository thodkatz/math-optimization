function a = bisection_curvature(curvature_condition, f, f_grad, x, pk, a, c=[1e-4, 0.9], ck, iter_count=20)
    % Implementation of bisection curvature condition (wolfe and goldstein)
    %
    % Note:
    % The algorithm is meant for wolfe weak. But analoguous, since goldstein and wolfe specify
    % a curvature condition, we use the same algorithm for goldstein
    %
    % Convergence of descent methods with backtracking 
    % (Armijo) linesearch. Bisection algorithm for weak Wolfe
    % conditions, Anton Evgrafov

    c1 = c(1);
    c2 = c(2);

    amin = 0;
    amax = inf;

    is_goldstein = strcmp(curvature_condition,"goldstein");
    is_wolfe = strcmp(curvature_condition,"weak wolfe");
    if not(is_goldstein || is_wolfe)
        error(-1)
    end

    iter = 0;
    % fprintf("\nSTARTED WOLFE Step size config...\n")
    phi_derivative0 = dot(f_grad(x), pk);
    fval = f(x);
    while 1
        if is_goldstein
            curv_cond = curvature_condition_goldstein(f, f_grad, fval, x, pk, a, c2, phi_derivative0);
        elseif is_wolfe
            curv_cond = curvature_condition_weak_wolfe(f_grad, x, pk, a, c2, phi_derivative0);
        else
            error(-1)
        end

        if armijo_condition(f, phi_derivative0, x, pk, a, c1, ck) 
            amax = a;
            a = (amin + amax)/2;
        elseif curv_cond
            amin = a;
            if amax == Inf;
                a *= 2;
            else
                a = (amin+amax)/2;
            end
        else
            break
        end
        iter += 1;
        if iter == iter_count
            fprintf("Bisection wolfe weak max iterations count. Couldn't satisfy conditions\n")
            break
        end
    end
    % fprintf("ENDED WOLFE Step size a=%f \n\n", a)
end

function b = curvature_condition_weak_wolfe(f_grad, x, pk, a, c, phi0)
    b = dot(f_grad(x + a*pk), pk) < c * phi0;
end

function b = curvature_condition_goldstein(f, f_grad, fval, x, pk, a, c, phi0)
    b = f(x + a*pk) < fval + (1-c)*a*phi0;
end
