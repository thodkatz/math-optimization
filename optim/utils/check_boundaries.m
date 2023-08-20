function a = check_boundaries(x, a, pk, domains)
    check_domain = x + a*pk;
    % clip out of bounds dimensions
    alpha_conditions = {};
    syms alpha
    to_check = false;
    for idx=1:numel(check_domain)
        bounds = domains{idx};
        lower_bound = bounds(1);
        upper_bound = bounds(2);
        xith = x(idx);
        pkith = pk(idx);
        check_domain_ith = check_domain(idx);
        if check_domain_ith < lower_bound
            if pkith < 0
                alpha_conditions(end+1) = alpha < (lower_bound - xith) / pkith;
                to_check = true;
            else
                alpha_conditions(end+1) = alpha > (lower_bound - xith) / pkith;
                to_check = true;
            end
        elseif check_domain_ith > upper_bound
            if pkith < 0
                alpha_conditions(end+1) = alpha > (upper_bound - xith) / pkith;
                to_check = true;
            else
                alpha_conditions(end+1) = alpha < (upper_bound - xith) / pkith;
                to_check = true;
            end
        end
    end
    if not(to_check)
        return
    end
    alpha_conditions(end+1) = alpha > 0;
    ret = solve(alpha_conditions{:});

    # parse lb and ub
    info = functions(function_handle(ret));
    splt = strsplit(info.function, "<");
    assert(numel(splt) == 3)
    lb = strrep(splt{1,1}, '@(alpha)', '');
    lb = strrep(lb, ' ', '');
    ub = strrep(strtrim(splt{1,3}), ' ', '');
    lb = str2num(lb);
    ub = str2num(ub);
    if strcmp(ub, 'inf')
        a = 1;
    else
        a = (ub - lb)/2;
    end
    fprintf("Checking boundaries... Alpa should be trimmed, new value %0.4e\n",a)

end
