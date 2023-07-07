function [a, num_fun_evals, num_grad_fun_evals] = bisection_curvature(curvature_condition, f, f_grad, x, pk, options, ck, num_fun_evals, num_grad_fun_evals)
    % Implementation of bisection curvature condition (wolfe and goldstein)
    %
    % Note:
    % The algorithm is meant for wolfe weak. But analoguous, since goldstein and wolfe specify
    % a curvature condition, we use the same algorithm for goldstein
    %
    % Convergence of descent methods with backtracking 
    % (Armijo) linesearch. Bisection algorithm for weak Wolfe
    % conditions, Anton Evgrafov

    a = options.a;
    iter_count = options.max_iters_step_size;

    is_goldstein = strcmp(curvature_condition,"goldstein");
    is_wolfe = strcmp(curvature_condition,"weak wolfe");

    if is_goldstein
        c1 = options.c;
        c2 = 1-c1;
        assert (c1>0 && c1<0.5)
    elseif is_wolfe
        c1 = options.c1;
        c2 = options.c2;
    else
        error(-1)
    end

    amin = 0;
    amax = inf;

    if not(is_goldstein || is_wolfe)
        error(-1)
    end

    iter = 0;
    fprintf("\nSTARTED WOLFE Step size config...\n")
    [fgrad, num_grad_fun_evals] = f_grad(x, num_grad_fun_evals);
    phi_derivative0 = dot(fgrad, pk);
    [fval, num_fun_evals]  = f(x, num_fun_evals);
    while 1
        if is_goldstein
            [curv_cond, num_fun_evals] = curvature_condition_goldstein(f, f_grad, fval, x, pk, a, c2, phi_derivative0, num_fun_evals);
        elseif is_wolfe
            [curv_cond, num_grad_fun_evals] = curvature_condition_weak_wolfe(f_grad, x, pk, a, c2, phi_derivative0, num_grad_fun_evals);
        else
            error(-1)
        end

        [armijo_cond, num_fun_evals] = armijo_condition(f, phi_derivative0, x, pk, a, c1, ck, num_fun_evals);
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
        fprintf("Iter %d, step size: %0.4e\n", iter, a)
        iter += 1;
        if iter == iter_count
            fprintf("Bisection wolfe weak max iterations count. Couldn't satisfy conditions\n")
            break
        end
    end
    fprintf("ENDED WOLFE Step size a=%0.4e \n\n", a)
end

function [b, number_of_evaluations] = curvature_condition_weak_wolfe(f_grad, x, pk, a, c, phi0, number_of_evaluations)
    [fgrad, number_of_evaluations] = f_grad(x+a*pk, number_of_evaluations);
    b = dot(fgrad, pk) < c * phi0;
end

function [b, number_of_evaluations] = curvature_condition_goldstein(f, f_grad, fval, x, pk, a, c, phi0, number_of_evaluations)
    [phi, number_of_evaluations] = f(x+a*pk, number_of_evaluations);
    b = phi < fval + c*a*phi0;
end
