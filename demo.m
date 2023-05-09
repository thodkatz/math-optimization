format longE

% setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

% define rosenbrock function
function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end


%search_x = -1.2:0.1:1.2;
%search_y = -1.2:0.1:1.2;
%[xmin, fmin] = newton(rosen_sym, [-1.8,-1.8]', "nonmonotone_armijo", search_x, search_y);

%function f = f2_sym()
    %syms x y
    %f = (x+2*y-7)^2 + (2*x + y - 5)^2;
%end

% search_x = 0:0.2:3;
% search_y = 0:0.2:3;
% [xmin, fmin] = steepest_descent(f2_sym, [3,3]', "none", search_x, search_y);

% search_x = -1.2:0.1:1.2;
% search_y = -1.2:0.1:1.2;
% c = [1e-4 0.9];
% rho = 2;
% [xmin, fmin] = steepest_descent(rosen_sym, [1.2,1.2]', "wolfe_strong", search_x, search_y, c, rho);

%function f = f4_sym()
    %syms x y
    %f = (x^2)^(y^2 + 1) + (y^2)^(x^2 + 1);
%end

% search_x = -1.5:0.5:1.5;
% search_y = -1.5:0.5:1.5;
% [xmin, fmin] = newton(f4_sym, [-1.5,1.25]', "backtracking_wolfe_weak", search_x, search_y);

% search_x = -1.5:0.5:1.5;
% search_y = -1.5:0.5:1.5;
% [xmin, fmin] = steepest_descent(f4_sym, [-1.5,1.25]', "backtracking_wolfe_weak", search_x, search_y);


%function f = f5_sym()
    %syms x y z
    %f = (x^2 + y^3 - z^4)^2 + (2*x*y*z)^2 + (2*x*y-3*y*z+x*z)^2;
%end

%search_x = -1:1:1;
%search_y = -1:1:1;
%[xmin, fmin] = newton(f5_sym, [10,10,10]', "nonmonotone_backtracking_wolfe_weak", search_x, search_y);

function f = f1_sym()
    syms x y
    f = x^2 + 4*y^2 + 2*x*y;
end

search_x = -3:0.2:0.4;
search_y = -3:0.2:0.4;
memory_limit = 100;
[xmin, fmin] = newton(f1_sym, [-3,-3]', 'grippo_bisection_wolfe_weak', search_x, search_y, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100, to_plot=false,memory_limit);
[xmin, fmin, iter] = steepest_descent(f1_sym, [-3,-3]', 'grippo_bisection_wolfe_weak', search_x, search_y, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100, to_plot=false,memory_limit);