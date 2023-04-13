function val = multidim_hessian(f_sym, point)
    % MULTIDIM_HESSIAN - Evaluate the Hessian matrix of a symbolic function at a point
    %
    % USAGE:
    %   val = multidim_hessian(f_sym, point)
    %
    % INPUTS:
    %   f_sym - symbolic function to evaluate
    %   point - point at which to evaluate the Hessian matrix
    %
    % OUTPUTS:
    %   val - Hessian matrix of f_sym evaluated at point
    %
    % DESCRIPTION:
    %   This function evaluates the Hessian matrix of a given symbolic function f_sym
    %   at a given point, using the built-in OCTAVE function `hessian`. The Hessian
    %   matrix is a matrix of second-order partial derivatives of f_sym with respect
    %   to its variables.
    %
    %   NOTE: This function requires the Symbolic Math Toolbox.
    %
    % EXAMPLES:
    %
    %   % Evaluate the Hessian of a quadratic function at a point
    %   syms x y
    %   f_sym = x^2 + 3*x*y + y^2;
    %   point = [1, 2];
    %   hessian_matrix = multidim_hessian(f_sym, point);

    [sym_vars, sym_point] = point2sym(f_sym, point);
    val = eval(subs(hessian(f_sym, sym_vars), sym_vars, sym_point));
end
