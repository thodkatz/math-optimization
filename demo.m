format longE

% setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

% define rosenbrock function
function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end


is_backtracking_wolfe_weak = 1;
is_wolfe_strong = 2;
is_backtracking_armijo = 3;
is_none = 4;

[xmin, fmin] = newton(rosen_sym, [1.2,1.2]', is_backtracking_wolfe_weak);

