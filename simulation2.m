%% Simulation for the functions listed on Chapter 4 of the doc

format short e

# setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

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

function [f,startp,xmin,fmin,domains] = rosen10()
    [f,startp,xmin,fmin,domains] = rosen_extended_sym(10);
end

function [f,startp,xmin,fmin,domains] = rosen20()
    [f,startp,xmin,fmin,domains] = rosen_extended_sym(20);
end

function [f,startp,xmin,fmin,domains] = wood()
    syms x1 x2 x3 x4 real
    f = 100*(x1^2-x2)^2 + (x1-1)^2 + (x3-1)^2 + 90*(x3^2 - x4)^2 + 10.1*((x2-1)^2 + (x4-1)^2) + 19.8*(x2-1)*(x4-1);
    startp = [-3,-1,-3,-1]';
    xmin = [1,1,1,1]';
    fmin = 0;
    domains = {[-inf,inf],[-inf,inf],[-inf,inf],[-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = powell()
    syms x1 x2 x3 x4 real
    f = (x1+10*x2)^2 + 5*(x3-x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^4;
    startp = [3,-1,0,1]';
    xmin = [0,0,0,0]';
    fmin = 0;
    domains = {[-inf,inf],[-inf,inf],[-inf,inf],[-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = cube()
    syms x1 x2 real
    f = 100*(x2-x1^3)^2 + (1-x1)^2
    startp = [-1.2,1]';
    xmin = [1,1]';
    fmin = 0;
    domains = {[-inf,inf],[-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = trig(n)
    syms x1
    f = n + (1-cos(x1)) - sin(x1) - cos(x1);
    startp = [1/(5*n)];
    xmin = [0];
    sym_vars = {x1};
    domains = {[-inf,inf]};

    for j=2:n
        sym_var = sym(['x', num2str(j)]);
        sym_vars{end+1} = sym_var;
        f -= cos(sym_var);
        startp(end+1) = 1/(5*n);
        xmin(end+1) = 0;
        domains{end+1} = [-inf,inf];
    end

    f = f^2;

    for i=2:n
        sym_var = sym_vars{i};
        ftemp = n + i*(1-cos(sym_var)) - sin(sym_var) - cos(sym_var);
        for j=2:n
            sym_var = sym_vars{j};
            ftemp -= cos(sym_var);
        end
        f += ftemp^2;
    end

    startp = startp';
    xmin = xmin';
    fmin = 0;
end

function [f,startp,xmin,fmin,domains] = trig20()
    [f,startp,xmin,fmin,domains] = trig(20);
end

function [f,startp,xmin,fmin,domains] = trig60()
    [f,startp,xmin,fmin,domains] = trig(60);
end


function [f1,f2,startp,xmin,fmin,domains] = helical()
    syms x1 x2 x3
    theta = atan(x2/x1)/(2*pi);
    f1 = 100*((x3-10*theta)^2 + (sqrt(x1^2+x2^2)-1)^2) + x3^2;

    theta += 1/2;
    f2 = 100*((x3-10*theta)^2 + (sqrt(x1^2+x2^2)-1)^2) + x3^2;
    startp = [-1,0,0]';
    xmin = [1,0,0]';
    fmin = 0;
    domains = {[-inf,inf],[-inf,inf],[-inf,inf]};
end



functions_info = {@rosen2,
            @rosen10,
            @rosen20,
            @wood, 
            @powell, 
            @cube, 
            @trig20};



line_search_methods = { 'backtracking_armijo', 
                        'wolfe_strong', 
                        'bisection_wolfe_weak',
                        'bisection_goldstein',
                        'hanger_zhang_backtracking_armijo', 
                        'hanger_zhang_bisection_wolfe_weak',
                        'hanger_zhang_wolfe_strong',
                        'hanger_zhang_bisection_goldstein',
                        'grippo_backtracking_armijo',
                        'grippo_bisection_wolfe_weak',
                        'grippo_wolfe_strong',
                        'grippo_bisection_goldstein'};



table = simulation(functions_info, line_search_methods);

dfile = 'simulation2.csv';
if exist(dfile, 'file') ; delete(dfile); end
cell2csv(dfile,table)


