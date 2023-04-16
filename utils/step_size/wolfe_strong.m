function a = wolfe_strong(f, f_grad, x, pk, a_max, rho = 2, c1=1e-4, c2=0.9, max_iters=20)
    a0 = 0
    a = a_max / 2
    a_prev = a0

    phi0 = f(x)
    phi_derivative0 = f_grad(x)
    phi_prev = phi0

    iter = 0
    fprintf("STARTED Wolfe strong line search...\n")
    while iter < max_iters
        phi = f(x + a*pk)
        if phi > phi0 + c1*a*phi_derivative0 || (iter > 0 && phi > phi_prev)
            a = zoom(a_prev, a)
            break
        end
        phi_derivative = f_grad(x+a*pk)
        if norm(phi_derivative) <= -c2*phi_derivative0
            break
        elseif phi_derivative >= 0
            a = zoom(a_prev, a)
            break
        end

        a_prev = a
        a = rho * a
        if a > a_max:
            a = a_max
        end

        phi_prev = phi

        fprintf("Iter %d, step size: %d\n", iter, a)
        iter += 1
    end

    if iter == max_iters:
        fprintf("Maximum number of iterations\n")
    end

    fprintf("ENDED Wolfe strong line search...\n")
end

function a_star = zoom(a_low, a_hi, phi_low, phi_high, phi0, derivative_phi0, x, f, f_grad, pk, c1=1e-4, c2=0.9, max_iters=10)
    % implementation: https://github.com/scipy/scipy/blob/v1.6.3/scipy/optimize/linesearch.py#L193-L335
    assert(a_low < a_hi)

    iter = 0
    delta1 = 0.2
    delta2 = 0.1
    phi_rec = phi0
    a_rec = 0
    fprintf("STARTED Zoom search")
    derivative_phi_low = f_grad(x+a_low*pk)
    while iter < max_iters
        dalpha = a_hi - a_low

        if iter > 0
            cubic_check = delta1 * dalpha
            a_j = cubicmin(a_low, phi_low, derivative_phi_low, a_hi, phi_high, a_rec, phi_rec)
        end

        if iter == 0 || a_j == null || a_hi - cubic_check < a_j  || a_j < a_low + cubic_check
            quad_check = delta2 * dalpha
            a_j = quadmin(a_low, phi_low, derivative_phi_low, a_hi, phi_high)
            if a_j == null || a_hi - quad_check < a_j || a_j < a_low + quad_check
                a_j = a_low + 0.5*dalpha
            end
        end

        phi_j = f(x+a_j*pk)
        if phi_j > phi0 + c1*a_j*derivative_phi0 || phi_j >= phi_low
            phi_rec = phi_high
            a_rec = a_hi
            a_hi = a_j
            phi_high = phi_j
        else
            derivative_phi_j = f_grad(x+a_j*pk)
            if norm(derivative_phi_j) <=  -c2*derivative_phi0
                a_star = a_j
                break
            elseif derivative_phi_j * (dalpha) >= 0
                phi_rec = phi_high
                a_rec = a_hi
                a_hi = a_low
                phi_high = phi_low
            else
                phi_rec = phi_low
                a_rec = a_low
            end
            a_low = a_j
            phi_low = phi_j
            derivative_phi_low = derivative_phi_j
        end
        fprintf("Iter %d, step size: %d\n", iter, a_j)
        iter += 1
    end

    if iter == max_iters
        fprintf("Maximum number of iterations\n")
        a_star = null
    end
    fprintf("ENDED Zoom search")
end

function cubicmin(x1, x2, y1, y2, x)

end

function quadmin()

end