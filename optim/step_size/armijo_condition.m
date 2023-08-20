function [b, number_of_evaluations] = armijo_condition(f, phi_derivative0, x, pk, a, c, ck, number_of_evaluations)
    [phi, number_of_evaluations] = f(x+a*pk, number_of_evaluations);
    b = phi > ck + c*a*phi_derivative0;
end
