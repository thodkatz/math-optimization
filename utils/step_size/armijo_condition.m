function b = armijo_condition(f, phi_derivative0, x, pk, a, c, ck)
    b = f(x + a*pk) > ck + c*a*phi_derivative0;
end
