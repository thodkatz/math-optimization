function val = sym2fun(f_sym, v)
    [sym_var, sym_point] = point2sym(f_sym, v);
    val = eval(subs(f_sym, sym_var, sym_point));
end