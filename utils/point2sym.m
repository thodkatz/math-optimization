function [sym_vars, sym_point] = point2sym(f_sym, point)
    % Convert a vector to a symbolic vector of floats to avoid the warning when substituting a symbolic function:
    % 'warning: passing floating-point values to sym is dangerous, see "help sym"'
    assert(is_vector(point))
    sym_vars = symvar(f_sym);
    assert(numel(sym_vars) == numel(point))

    sym_point = [];
    for i = 1:numel(sym_vars)
        sym_point = [sym_point, sym(point(i), "f")];
    end
end