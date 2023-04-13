function val = multidim_grad(f_sym, point)
    % MULTIDIM_GRAD - Evaluate the gradient of a symbolic function at a point
    %
    % USAGE:
    %   val = multidim_grad(f_sym, point)
    %
    % INPUTS:
    %   f_sym - symbolic function to evaluate
    %   point - point at which to evaluate the gradient
    %
    % OUTPUTS:
    %   val - gradient of f_sym evaluated at point
    %
    % DESCRIPTION:
    %   This function evaluates the gradient of a given symbolic function f_sym at a
    %   given point, using the built-in OCTAVE function `gradient`. The gradient is a
    %   vector of first-order partial derivatives of f_sym with respect to its
    %   variables.
    %
    %   NOTE: This function requires the Symbolic Math Toolbox.
    %
    % EXAMPLES:
    %
    %   % Evaluate the gradient of a quadratic function at a point
    %   syms x y
    %   f_sym = x^2 + 3*x*y + y^2;
    %   point = [1, 2];
    %   gradient_vector = multidim_grad(f_sym, point);

    [sym_vars, sym_point] = point2sym(f_sym, point);
    val = eval(subs(gradient(f_sym, sym_vars), sym_vars, sym_point));
end