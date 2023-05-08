%% Simulation for the functions listed on Chapter 3 of the doc

format longE

# setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

function [e] = get_error(expected, observed)
    % if the observed value isn't real, then the error should be Inf
    if isreal(observed)
        e = abs(expected - observed);
    else
        e = Inf;
    end
end

function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end

function f = f1_sym()
    syms x y
    f = x^2 + 4*y^2 + 2*x*y;
end

function f = f2_sym()
    syms x y
    f = (x+2*y-7)^2 + (2*x + y - 5)^2;
end

function f = f3_sym()
    syms x y 
    f = 5*x^4 + 6*y^4 - 6*x^2 + 2*x*y + 5*y^2 + 15*x - 7*y + 13;
end

function f = f4_sym()
    syms x y
    f = (x^2)^(y^2 + 1) + (y^2)^(x^2 + 1);
end

function f = f5_sym()
    syms x y z
    f = (x^2 + y^3 - z^4)^2 + (2*x*y*z)^2 + (2*x*y-3*y*z+x*z)^2;
end


function f = f6_sym()
    syms x y z k
    f = (x-1)^2 + (x-sqrt(y))^2 + (y-sqrt(z))^2 + (z-sqrt(k))^2;
end


functions = {rosen_sym, 
            f1_sym, 
            f2_sym, 
            f3_sym, 
            f4_sym, 
            f5_sym, 
            f6_sym};

starting_points = {[-1.8, -1.8]',
                    [-3,-3]',
                    [-9.5,9.5]',
                    [1.9,-1.9]',
                    [-1.5,1.25]',
                    [10,10,10]',
                    [0.1,0.1,0.1,0.1]'}
expected_fmin = {0,
                0,
                0,
                -6.4931,
                0,
                0,
                1.13719e-10};
expected_xmin = {[1,1]',
                    [0,0]',
                    [1,3]',
                    [-1.1515,0.5455]',
                    [0,0]',
                    [0.0048605831,0.0016994507, 0]',
                    [0.999993, 0.999983, 0.999964, 0.999912]'};

to_plot=false;
% redundant search_x and search_y because we don't plot
% todo use nargin, to configure the inputs of the functions better
search_x = -1.8:0.1:1.2;
search_y = -1.8:0.1:1.2;

line_search_methods = {'none', 
                        'backtracking_armijo', 
                        'wolfe_strong', 
                        'bisection_wolfe_weak',
                        'nonmonotone_backtracking_armijo', 
                        'nonmonotone_bisection_wolfe_weak',
                        'nonmonotone_wolfe_strong'};

 table = {'function', 'method', 'step size', 'iterations', 'error x1', 'error x2', 'error fvalue'};

count = 1
for i=1:numel(functions)
    fun = functions{i}
    starting_point = starting_points{i}
    for j=1:numel(line_search_methods)
        line_search_method = line_search_methods{j}
        if strcmp(line_search_method, 'wolfe_strong') || strcmp(line_search_method, 'nonmonotone_wolfe_strong')
            c = [1e-4 0.9];
            rho = 2;
        else
            c = 0.1;
            rho = 0.5;
        end
        [xmin, fmin, iter] = newton(fun, 
                            starting_point, 
                            line_search_method, 
                            search_x, 
                            search_y, 
                            c, 
                            rho, 
                            a=1, 
                            eps=1e-6, 
                            max_iters=100, 
                            to_plot=to_plot);

        expected_x = expected_xmin{i};
        table(end+1,:) = {i, "newton", line_search_method, iter, get_error(expected_x(1),xmin(1)), get_error(expected_x(2),xmin(2)), get_error(expected_fmin{i},fmin)};
    end
    for j=1:numel(line_search_methods)
        try
            [xmin, fmin, iter] = steepest_descent(fun, 
                                starting_point, 
                                line_search_method, 
                                search_x, 
                                search_y, 
                                c, 
                                rho, 
                                a=1, 
                                eps=1e-6, 
                                max_iters=100, 
                                to_plot=to_plot);
            expected_x = expected_xmin{i};
            table(end+1,:) = {i, "steepest descent", line_search_method, iter, get_error(expected_x(1),xmin(1)), get_error(expected_x(2),xmin(2)), get_error(expected_fmin{i},fmin)};
        catch err
            err.message
            fprintf("Failed to converge\n")
            table(end+1,:) = {i, "steepest descent", line_search_method, iter, Inf, Inf, Inf};
        end_try_catch
    end
end

dataframe(table)

