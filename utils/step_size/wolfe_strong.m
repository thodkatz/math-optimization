function [a, number_of_function_evaluations, number_of_gradient_function_evaluations] = wolfe_strong(f, f_grad, x, pk, a_max, rho = 2, c=[1e-4, 0.9], ck, max_iters=20, number_of_function_evaluations, number_of_gradient_function_evaluations)
    % Jorge Nocedal and Stephen J Wright. Numerical optimization. Springer, 1999
    % pages: 60-61
    %
    % https://github.com/scipy/scipy/blob/v1.6.3/scipy/optimize/linesearch.py#L193-L335

    c1 = c(1);
    c2 = c(2);

    a0 = 0;
    a = a_max/2;
    a_prev = a0;

    [fgrad, number_of_gradient_function_evaluations] = f_grad(x, number_of_gradient_function_evaluations);
    phi_derivative0 = dot(fgrad,pk);
    phi_prev = ck;

    iter = 0;
    % fprintf("\nSTARTED Wolfe strong line search...\n")
    while iter < max_iters
        [phi, number_of_function_evaluations] = f(x + a*pk, number_of_function_evaluations);
        if phi > ck + c1*a*phi_derivative0 || (iter > 0 && phi > phi_prev)
            [a, number_of_function_evaluations, number_of_gradient_function_evaluations] = zoom(a_prev, a, phi_prev, phi, ck, phi_derivative0, x, f, f_grad, pk, c1, c2, 10, number_of_function_evaluations, number_of_gradient_function_evaluations);
            break
        end
        [fgrad, number_of_gradient_function_evaluations] = f_grad(x+a*pk);
        phi_derivative = dot(fgrad,pk);
        if norm(phi_derivative) <= -c2*phi_derivative0
            break
        elseif phi_derivative >= 0
            [a, number_of_function_evaluations, number_of_gradient_function_evaluations] = zoom(a_prev, a, phi_prev, phi, ck, phi_derivative0, x, f, f_grad, pk, c1, c2, 10, number_of_function_evaluations, number_of_gradient_function_evaluations);
            break
        end

        a_prev = a;
        a = rho * a;
        if a > a_max
            a = a_max;
        end

        phi_prev = phi;

        % fprintf("Iter %d, step size: %d\n", iter, a)
        iter += 1;
    end

    if iter == max_iters
        % fprintf("Maximum number of iterations\n")
    end

    % fprintf("ENDED Wolfe strong line search a=%f \n\n", a)
end

function [a_star, number_of_function_evaluations, number_of_gradient_function_evaluations] = zoom(a_low, a_hi, phi_low, phi_high, ck, derivative_phi0, x, f, f_grad, pk, c1=1e-4, c2=0.9, max_iters=10, number_of_function_evaluations, number_of_gradient_function_evaluations)
    % implementation: https://github.com/scipy/scipy/blob/v1.6.3/scipy/optimize/linesearch.py#L193-L335
    assert(a_low < a_hi)

    iter = 0;
    delta1 = 0.2;
    delta2 = 0.1;
    % fprintf("STARTED Zoom search\n")
    [derivative_phi_low, number_of_gradient_function_evaluations] = f_grad(x+a_low*pk, number_of_gradient_function_evaluations);
    while iter < max_iters
        dalpha = a_hi - a_low;

        # bisection interpolation
        a_j = a_low + 0.5*dalpha;

        [phi_j, number_of_function_evaluations] = f(x+a_j*pk, number_of_function_evaluations);
        if phi_j > ck + c1*a_j*derivative_phi0 || phi_j >= phi_low
            a_hi = a_j;
            phi_high = phi_j;
        else
            [derivative_phi_j, number_of_gradient_function_evaluations] = f_grad(x+a_j*pk, number_of_gradient_function_evaluations);
            if norm(derivative_phi_j) <=  -c2*derivative_phi0
                a_star = a_j;
                break
            elseif derivative_phi_j * (dalpha) >= 0
                a_hi = a_low;
                phi_high = phi_low;
            end
            a_low = a_j;
            phi_low = phi_j;
            derivative_phi_low = derivative_phi_j;
        end
        % fprintf("Iter %d, step size: %d\n", iter, a_j)
        iter += 1;
    end

    if iter == max_iters
        % fprintf("Maximum number of iterations\n")
        a_star = a_j; % not sure what to do here?
    end
    % fprintf("ENDED Zoom search a_star = %f \n", a_star)
end
