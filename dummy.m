format longE

% setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

% define rosenbrock function
% function f = rosen_sym()
%     syms x y
%     f = 100*(y-x^2)^2 + (1-x)^2;
% end


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

function f = f7_sym()
    syms x1 x2 x3
    f = x1^2 + (x2 + x2^2)^2 + (-1 + exp(x3))^2;
end

function [f,domains] = f8_sym()
    syms x1 x2 x3 x4
    f = (x1-1)^2 + (x1-sqrt(x2))^2 + (x2-sqrt(x3))^2 + (x3-sqrt(x4))^2;
    domains = {[-inf,inf], [0,inf], [0,inf], [0,inf]}; # exclude zero take into account grad
end

% search_x = -1:1:1;
% search_y = -1:1:1;
% c = [1e-4 0.9];
% rho = 2;
% [xmin, fmin] = newton(f7_sym, [2,3,-8]', 'wolfe_strong', search_x, search_y, c, rho);

function f = f1_sym()
    syms x y
    f = x^2 + 4*y^2 + 2*x*y;
end

% search_x = -3:0.2:0.4;
% search_y = -3:0.2:0.4;
% memory_limit = 100;

% c = [1e-4 0.9];
% rho = 2;
% [xmin, fmin] = steepest_descent(rosen_sym, [-1.8,-1.8]', 'hanger_zhang_wolfe_strong', search_x, search_y, c=c, rho=rho, a=1, eps=1e-6, max_iters=1000, to_plot=false);


% [xmin, fmin] = newton(f1_sym, [-3,-3]', 'grippo_bisection_wolfe_weak', search_x, search_y, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100, to_plot=false,memory_limit);
% [xmin, fmin, iter] = steepest_descent(f1_sym, [-3,-3]', 'grippo_bisection_wolfe_weak', search_x, search_y, c=0.1, rho=0.5, a=1, eps=1e-6, max_iters=100, to_plot=false,memory_limit);

search_x = -1:1:1;
search_y = -1:1:1;
c = [1e-4, 0.9];
rho = 2;
a = 1;
eps = 1e-6;
max_iters=4;
[f8,domains] = f8_sym
% [xmin, fmin] = newton(f8, domains, [0.1,0.1,0.1,0.1]', 'bisection_wolfe_weak', search_x, search_y, c=c, rho=rho,a=a,eps,max_iters=max_iters, to_plot=false);

function [f,startp,xmin,fmin,domains] = rosen_extended_sym(n)
    x1 = sym('x1');
    x2 = sym('x2');
    f = 100*(x2-x1^2)^2 + (1-x1)^2;

    startp = [-1.2,1];
    xmin = [1,1];
    domains = {[-inf,inf], [-inf,inf]};

    for i = 3:n
        sym_ith = sym(['x',num2str(i)]);
        if i == 3
            sym_ith_prev = x2;
        else
            sym_ith_prev = sym(['x',num2str(i-1)]);
        end
        f += (1 - sym_ith_prev)^2 + 100*(sym_ith - sym_ith_prev^2)^2;

        if mod(i,2) == 0
            startp(end+1) = 1;
        else
            startp(end+1) = -1.2;
        end
        xmin(end+1) = 1;
        domains{end+1} = [-inf,inf];
    end
    startp = startp';
    xmin = xmin';
    fmin = 0;
end

function [f,startp,xmin,fmin,domains] = rosen2()
    [f,startp,xmin,fmin,domains] = rosen_extended_sym(2);
end

function [f,startp,xmin,fmin,domains] = wood()
    syms x1 x2 x3 x4 real
    f = 100*(x1^2-x2)^2 + (x1-1)^2 + (x3-1)^2 + 90*(x3^2 - x4)^2 + 10.1*((x2-1)^2 + (x4-1)^2) + 19.8*(x2-1)*(x4-1);
    startp = [-3,-1,-3,-1]';
    xmin = [1,1,1,1]';
    fmin = 0;
    domains = {[-inf,inf],[-inf,inf],[-inf,inf],[-inf,inf]};
end

config.a = 1;
config.rho = 0.6;
config.c = 0.1;
config.eps=1e-16;
config.max_iters = 100;
config.max_iters_step_size = 50;
config.memory_limit = 0;
[f,startp,xmin,fmin,domains] = wood();
% [xmin, fmin] = newton(f, domains, startp,"grippo_backtracking_armijo",config);

function [f,startp,xmin,fmin,domains] = powell()
    syms x1 x2 x3 x4 real
    f = (x1+10*x2)^2 + 5*(x3-x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4;
    startp = [3,-1,0,1]';
    xmin = [0,0,0,0]';
    fmin = 0;
    domains = {[-inf,inf],[-inf,inf],[-inf,inf],[-inf,inf]};
end
config.a = 1;
config.rho = 0.6;
config.c = 0.1;
config.eps=1e-16;
config.max_iters = 100;
config.max_iters_step_size = 50;
config.memory_limit = 0;
[f,startp,xmin,fmin,domains] = powell();
[xmin, fmin] = newton(f, domains, startp,"grippo_backtracking_armijo",config);

