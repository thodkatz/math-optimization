function table = simulation(functions_info, line_search_methods)
    table = {'function', 'method', 'step size', 'iterations', 'ni', 'nj', 'error x1', 'error x2', 'error fvalue'};
    to_plot=false;
    % redundant search_x and search_y because we don't plot
    % todo use nargin, to configure the inputs of the functions better
    search_x = -1.8:0.1:1.2;
    search_y = -1.8:0.1:1.2;
    config.a=1;
    config.eps = 1e-16;
    config.memory_limit = 10; % for grippo
    config.max_iters = 5e3;
    config.max_iters_step_size = 50;
    config.max_iters_zoom = 10; % for wolfe strong

    for i=1:numel(functions_info)
        fun_info = functions_info{i};
        [fsym, starting_point, expected_x, expected_f,domains] = fun_info()
        for j=1:numel(line_search_methods)
            line_search_method = line_search_methods{j}

            if endsWith(line_search_method, 'wolfe_strong')
                config.c1 = 1e-4;
                config.c2 = 0.9;
                config.rho = 2;
            else
                config.c = 0.1;
                config.rho = 0.6;
            end

            [xmin, fmin, iter, ni, nj] = newton(fsym, domains, starting_point, line_search_method, config);

            table(end+1,:) = {['f',num2str(i)], "newton", line_search_method, iter, ni, nj, get_error(expected_x(1),xmin(1)), get_error(expected_x(2),xmin(2)), get_error(expected_f,fmin)};
        end
        for j=1:numel(line_search_methods)
            line_search_method = line_search_methods{j}

            if endsWith(line_search_method, 'wolfe_strong')
                config.c1 = 1e-4;
                config.c2 = 0.9;
                config.rho = 2;
            else
                config.c = 0.1;
                config.rho = 0.6;
            end

            try
                [xmin, fmin, iter, ni, nj] = steepest_descent(fsym, domains, starting_point, line_search_method, config);

                table(end+1,:) = {['f',num2str(i)], "steepest", line_search_method, iter, ni, nj, get_error(expected_x(1),xmin(1)), get_error(expected_x(2),xmin(2)), get_error(expected_f,fmin)};
            catch err
                err.message
                fprintf("Failed to converge\n")
                table(end+1,:) = {['f',num2str(i)], "steepest", line_search_method, iter, ni, nj, Inf, Inf, Inf};
            end_try_catch
        end
    end

    % format table
    [rows, cols] = size(table);
    for i=2:rows
        for j=1:cols
        header = table(1,j);
            if strcmp(header, "iterations") || strcmp(header, "ni") || strcmp(header, "nj")
                table{i,j} = sprintf("%d", table{i,j});
            end
        end
    end
end

function [e] = get_error(expected, observed)
    % if the observed value isn't real, then the error should be Inf
    if isreal(observed)
        e = abs(expected - observed);
    else
        e = Inf;
        error(-1);
    end
end