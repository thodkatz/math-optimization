format longE

# setup the path to include the 'utils' directory
directory = fileparts(which(mfilename)); 
addpath(genpath(directory));    


function main
    is_wolfe = 1;
    is_armijo = 2;
    is_none = 3;


    % run_rosen(is_armijo)
    run_rosen(is_wolfe)
    run_rosen(is_none)
end

function run_rosen(line_searh_method)
    [xmin, ymin] = newton(@rosen, @rosen_grad, [1.2,1.2]', line_searh_method)

    [xmin, ymin] = newton(@rosen, @rosen_grad, [-1.2,1]', line_searh_method)

    [xmin, ymin] = steepest_descent(@rosen, @rosen_grad, [1.2,1.2]', line_searh_method)

    [xmin, ymin] = steepest_descent(@rosen, @rosen_grad, [-1.2,1]', line_searh_method)
end

function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end

function val = rosen(v)
    [sym_var, sym_point] = point2sym(rosen_sym, v);
    val = eval(subs(rosen_sym, sym_var, sym_point));
end

function val = rosen_grad(point)
    assert(is_vector(point))
    val = multidim_grad(rosen_sym, point);
end

function val = rosen_hessian(v)
    # Given a vector 'v' with 2 coordinates return the hessian of the Rosenbrock function, a 2x2 matrix.
    val = multidim_hessian(rosen_sym, v);
end


main()