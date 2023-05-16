function [a, number_of_function_evaluations, number_of_gradient_function_evaluations] = bisection_curvature(curvature_condition, f, f_grad, x, pk, a, c=[1e-4, 0.9], ck, iter_count=20, number_of_function_evaluations, number_of_gradient_function_evaluations)
    % Implementation of bisection curvature condition (wolfe and goldstein)
    %
    % Note:
    % The algorithm is meant for wolfe weak. But analoguous, since goldstein and wolfe specify
    % a curvature condition, we use the same algorithm for goldstein
    %
    % Convergence of descent methods with backtracking 
    % (Armijo) linesearch. Bisection algorithm for weak Wolfe
    % conditions, Anton Evgrafov
    is_goldstein = strcmp(curvature_condition,"goldstein");
    is_wolfe = strcmp(curvature_condition,"weak wolfe");

    if is_goldstein
        c1 = c;
        c2 = 1-c;
        assert (c1>0 && c1<0.5)
    elseif is_wolfe
        c1 = c(1);
        c2 = c(2);
    else
        error(-1)
    end

    amin = 0;
    amax = inf;

    if not(is_goldstein || is_wolfe)
        error(-1)
    end

    iter = 0;
    % fprintf("\nSTARTED WOLFE Step size config...\n")
    [fgrad, number_of_gradient_function_evaluations] = f_grad(x, number_of_gradient_function_evaluations);
    phi_derivative0 = dot(fgrad, pk);
    [fval, number_of_function_evaluations]  = f(x, number_of_function_evaluations);
    while 1
        if is_goldstein
            [curv_cond, number_of_function_evaluations] = curvature_condition_goldstein(f, f_grad, fval, x, pk, a, c2, phi_derivative0, number_of_function_evaluations);
        elseif is_wolfe
            [curv_cond, number_of_gradient_function_evaluations] = curvature_condition_weak_wolfe(f_grad, x, pk, a, c2, phi_derivative0, number_of_gradient_function_evaluations);
        else
            error(-1)
        end

        [armijo_cond, number_of_function_evaluations] = armijo_condition(f, phi_derivative0, x, pk, a, c1, ck, number_of_function_evaluations);
        if armijo_cond
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

function [b, number_of_evaluations] = curvature_condition_weak_wolfe(f_grad, x, pk, a, c, phi0, number_of_evaluations)
    [fgrad, number_of_evaluations] = f_grad(x+a*pk, number_of_evaluations);
    b = dot(fgrad, pk) < c * phi0;
end

function [b, number_of_evaluations] = curvature_condition_goldstein(f, f_grad, fval, x, pk, a, c, phi0, number_of_evaluations)
    [phi, number_of_evaluations] = f(x+a*pk, number_of_evaluations);
    b = phi < fval + c*a*phi0;
end
