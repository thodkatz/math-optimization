format longE

% setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

% define rosenbrock function
function [f,startp,xmin,fmin,domains] = f1_sym()
    syms x1 x2
    f = 100*(x2-x1^2)^2 + (1-x1)^2;
    startp = [-1.8,-1.8]';
    xmin = [1,1]';
    fmin = 0;
    domains = {[-inf,inf], [-inf,inf]};
end

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


config.a = 1;
config.rho = 0.5;
config.c = 0.1;
config.eps = 1e-16;
config.max_iters = 200;
config.max_iters_step_size = 100;
search_x = -1.2:0.1:1.2;
search_y = -1.2:0.1:1.2;
[f1,startp,xmin,fmin,domains] = f1_sym();

% fmin 0: (a=1,rho=0.5,c=0.1,eps=1e-16,max_iters=200,max_iters_step_size=100)
% [xmin, fmin] = newton(f1, domains, [-1.8,-1.8]',"backtracking_armijo",config);

% fmin 1e-6: (a=1,rho=0.5,c=0.1,eps=1e-16,max_iters=5000,max_iters_step_size=50), the last one doesn't matter, usually the backtracking finds a solution within 8 iters
% fmin 5.41e-8: (a=1,rho=0.6,c=0.1,eps=1e-16, max_iters=5000,max_iters_step_size=50)
% fmin 4.42e-12: (a=1,rho=0.6,c=0.1,eps=1e-16, max_iters=1e4,max_iters_step_size=50)
% fmin 2.25e-10: (a=1,rho=0.7,c=0.1,eps=1e-16,max_iters=1e4,max_iters_step_size=50)
% increasing the numbers of iterations it is converging
config.a = 1;
config.rho = 0.6;
config.c = 0.1;
config.eps=1e-16;
config.max_iters = 1000;
config.max_iters_step_size = 50;
% [xmin, fmin] = steepest_descent(f1, domains, [-1.8,-1.8]',"backtracking_armijo",config);

% fmin 3e-31 (a=1,rho=2,c1=1e-4,c2=0.9,eps=1e-16,max_iters=1000,max_iters_step_size=50,max_iters_zoom=10)
config.a = 1;
config.rho = 2;
config.c1 = 1e-4;
config.c2 = 0.9;
config.eps=1e-16;
config.max_iters = 1000;
config.max_iters_step_size = 50;
config.max_iters_zoom = 10;
% [xmin, fmin] = newton(f1, domains, [-1.8,-1.8]',"wolfe_strong",config);

% fmin 6.86e-03 (a=1,rho=2,c1=1e-4,c2=0.9,eps=1e-16,max_iters=1000,max_iters_step_size=50,max_iters_zoom=10)
% fmin 9.73e-04 (a=1,rho=2,c1=1e-4,c2=0.9,eps=1e-16,max_iters=1e4,max_iters_step_size=50,max_iters_zoom=10)
config.a = 1;
config.rho = 2;
config.c1 = 1e-4;
config.c2 = 0.9;
config.eps=1e-16;
config.max_iters = 1e4;
config.max_iters_step_size = 50;
config.max_iters_zoom = 10;
% [xmin, fmin] = steepest_descent(f1, domains, [-1.8,-1.8]',"wolfe_strong",config);

% fmin 0 (a=1,rho=2,c1=1e-4,c2=0.9,rho=2,eps=1e-16,max_iters=1e4,max_iters_step_size=50,max_iters_zoom=10,memory_limit=50)
config.a = 1;
config.rho = 2;
config.c1 = 1e-4;
config.c2 = 0.9;
config.eps=1e-16;
config.max_iters = 1e4;
config.max_iters_step_size = 50;
config.max_iters_zoom = 10;
config.memory_limit = 50;
% [xmin, fmin] = newton(f1, domains, [-1.8,-1.8]',"grippo_wolfe_strong");

config.a = 1;
config.rho = 5;
config.c1 = 1e-4;
config.c2 = 0.9;
config.eps=1e-16;
config.max_iters = 200;
config.max_iters_step_size = 50;
config.max_iters_zoom = 50;
config.memory_limit = 50;
% [xmin, fmin] = steepest_descent(f1, domains, [-1.8,-1.8]',"grippo_wolfe_strong",config);

% fmin 0 (a=1,rho=0.6,c=0.1,eps=1e-16,max_iters=1e4,max_iters_step_size=50,memory_limit=50)
config.a = 1;
config.rho = 0.6;
config.c = 0.1;
config.eps=1e-16;
config.max_iters = 1e4;
config.max_iters_step_size = 50;
config.memory_limit = 50;
% [xmin, fmin] = newton(f1, domains, [-1.8,-1.8]',"grippo_backtracking_armijo",config);

% fmin  4.76 (a=1,rho=0.6,c=0.1,eps=1e-16,max_iters=1e4,max_iters_step_size=50,memory_limit=50)
config.a = 1;
config.rho = 0.6;
config.c = 0.1;
config.eps=1e-16;
config.max_iters = 1000;
config.max_iters_step_size = 50;
config.memory_limit = 10;
[xmin, fmin] = steepest_descent(f1, domains, [-1.8,-1.8]',"grippo_backtracking_armijo",config,"search_x",search_x,"search_y",search_y);

%
config.a = 1;
config.rho = 2;
config.c1 = 1e-4;
config.c2 = 0.9;
config.eps=1e-16;
config.max_iters = 100;
config.max_iters_step_size = 50;
config.max_iters_zoom = 10;
config.memory_limit = 20;
% [xmin, fmin] = steepest_descent(rosen2, domains, [-1.2,1]',"grippo_wolfe_strong",config,"search_x",search_x,"search_y",search_y);
% [xmin, fmin] = steepest_descent(rosen2, domains, [-1.2,1]',"grippo_wolfe_strong",config,"search_x",search_x,"search_y",search_y);


% fmin 6.43 (a=1,rho=0.6,c=0.1,eps=1e-16,max_iters=1e4,max_iters_step_size=50)
config.a = 1;
config.rho = 0.6;
config.c = 0.1;
config.eps=1e-16;
config.max_iters = 1e4;
config.max_iters_step_size = 50;
% [xmin, fmin] = steepest_descent(f1, domains, [-1.8,-1.8]',"hanger_zhang_backtracking_armijo",config,"search_x",search_x,"search_y",search_y);