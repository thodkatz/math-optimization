function table = simulation(functions_info, line_search_methods)
    table = {'function', 'method', 'step size', 'iterations', 'ni', 'nj', 'error x1', 'error x2', 'error fvalue'};
    to_plot=false;
    % redundant search_x and search_y because we don't plot
    % todo use nargin, to configure the inputs of the functions better
    search_x = -1.8:0.1:1.2;
    search_y = -1.8:0.1:1.2;
    memory_limit = 50; % for grippo
    for i=1:numel(functions_info)
        fun_info = functions_info{i};
        [fsym, starting_point, expected_x, expected_f] = fun_info()
        for j=1:numel(line_search_methods)
            line_search_method = line_search_methods{j}
            if endsWith(line_search_method, 'wolfe_strong')
                c = [1e-4 0.9];
                rho = 2;
            else
                c = 0.1;
                rho = 0.5;
            end
            if startsWith(line_search_method, "grippo")
                [xmin, fmin, iter, ni, nj] = newton(fsym, starting_point, line_search_method, search_x, search_y, c, rho, a=1, eps=1e-6, max_iters=100, to_plot=to_plot, memory_limit);
            else
                [xmin, fmin, iter, ni, nj] = newton(fsym, starting_point, line_search_method, search_x, search_y, c, rho, a=1, eps=1e-6, max_iters=100, to_plot=to_plot);
            end

            table(end+1,:) = {['f',num2str(i)], "newton", line_search_method, iter, ni, nj, get_error(expected_x(1),xmin(1)), get_error(expected_x(2),xmin(2)), get_error(expected_f,fmin)};
        end
        for j=1:numel(line_search_methods)
            line_search_method = line_search_methods{j}
            try
                if startsWith(line_search_method, "grippo")
                    [xmin, fmin, iter, ni, nj] = steepest_descent(fsym, starting_point, line_search_method, search_x, search_y, c, rho, a=1, eps=1e-6, max_iters=100, to_plot=to_plot, memory_limit);
                else
                    [xmin, fmin, iter, ni, nj] = steepest_descent(fsym, starting_point, line_search_method, search_x, search_y, c, rho, a=1, eps=1e-6, max_iters=100, to_plot=to_plot);
                end

                table(end+1,:) = {['f',num2str(i)], "steepest", line_search_method, iter, ni, nj, get_error(expected_x(1),xmin(1)), get_error(expected_x(2),xmin(2)), get_error(expected_f,fmin)};
            catch err
                err.message
                fprintf("Failed to converge\n")
                table(end+1,:) = {['f',num2str(i)], "steepest", line_search_method, iter, ni, nj, Inf, Inf, Inf};
            end_try_catch
        end
    end

    % format table
    for i=1:numel(table)
        if isnumeric(table{i})
            table{i} = sprintf("%.4e", table{i});
        end
    end
end

function [e] = get_error(expected, observed)
    % if the observed value isn't real, then the error should be Inf
    if isreal(observed)
        e = abs(expected - observed);
    else
        e = Inf;
    end
end