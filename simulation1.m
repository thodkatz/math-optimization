%% Simulation for the functions listed on Chapter 3 of the doc

format short e

# setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))



function [f,startp,xmin,fmin] = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
    startp = [-1.8,-1.8]';
    xmin = [1,1]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = f1_sym()
    syms x y
    f = x^2 + 4*y^2 + 2*x*y;
    startp = [-3,-3]';
    xmin = [0,0]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = f2_sym()
    syms x y
    f = (x+2*y-7)^2 + (2*x + y - 5)^2;
    startp = [-9.5,-9.5]';
    xmin = [1,3]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = f3_sym()
    syms x y 
    f = 5*x^4 + 6*y^4 - 6*x^2 + 2*x*y + 5*y^2 + 15*x - 7*y + 13;
    startp = [1.9,-1.9]';
    xmin = [-1.1515,0.5455]';
    fmin = -6.4931;
end

function [f,startp,xmin,fmin] = f4_sym()
    syms x y
    f = (x^2)^(y^2 + 1) + (y^2)^(x^2 + 1);
    startp = [-1.5,11.25]';
    xmin = [0,0]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = f5_sym()
    syms x y z
    f = (x^2 + y^3 - z^4)^2 + (2*x*y*z)^2 + (2*x*y-3*y*z+x*z)^2;
    startp = [10,10,10]';
    xmin = [0.0048605831,0.0016994507, 0]';
    fmin = 0;
end


function [f,startp,xmin,fmin] = f6_sym()
    syms x y z k
    f = (x-1)^2 + (x-sqrt(y))^2 + (y-sqrt(z))^2 + (z-sqrt(k))^2;
    startp = [0.1,0.1,0.1,0.1]';
    xmin = [0.999993, 0.999983, 0.999964, 0.999912]';
    fmin = 1.13719e-10;
end


functions_info = {@rosen_sym, 
            @f1_sym, 
            @f2_sym, 
            @f3_sym, 
            @f4_sym, 
            @f5_sym, 
            @f6_sym};

starting_points = {[-1.8, -1.8]',
                    [-3,-3]',
                    [-9.5,9.5]',
                    [1.9,-1.9]',
                    [-1.5,1.25]',
                    [10,10,10]',
                    [0.1,0.1,0.1,0.1]'};
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



line_search_methods = {'none', 
                        'backtracking_armijo', 
                        'wolfe_strong', 
                        'bisection_wolfe_weak',
                        'bisection_goldstein',
                        'hanger_zhang_backtracking_armijo', 
                        'hanger_zhang_bisection_wolfe_weak',
                        'hanger_zhang_wolfe_strong',
                        'grippo_backtracking_armijo',
                        'grippo_bisection_wolfe_weak',
                        'grippo_wolfe_strong'};

table = simulation(functions_info, line_search_methods);

dfile = 'simulation1.csv';
if exist(dfile, 'file') ; delete(dfile); end
cell2csv(dfile,table)
