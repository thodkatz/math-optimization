%% Simulation for the functions listed on Chapter 3 of the doc

format short e

# setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))



function [f,startp,xmin,fmin,domains] = f1_sym()
    syms x1 x2
    f = 100*(x2-x1^2)^2 + (1-x1)^2;
    startp = [-1.8,-1.8]';
    xmin = [1,1]';
    fmin = 0;
    domains = {[-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f2_sym()
    syms x1 x2
    f = x1^2 + 4*x2^2 + 2*x1*x2;
    startp = [-3,-3]';
    xmin = [0,0]';
    fmin = 0;
    domains = {[-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f3_sym()
    syms x1 x2
    f = (x1+2*x2-7)^2 + (2*x1 + x2 - 5)^2;
    startp = [-9.5,-9.5]';
    xmin = [1,3]';
    fmin = 0;
    domains = {[-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f4_sym()
    syms x1 x2
    f = 5*x1^4 + 6*x2^4 - 6*x1^2 + 2*x1*x2 + 5*x2^2 + 15*x1 - 7*x2 + 13;
    startp = [1.9,-1.9]';
    xmin = [-1.1515,0.5455]';
    fmin = -6.4931;
    domains = {[-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f5_sym()
    syms x1 x2
    f = (x1^2)^(x2^2 + 1) + (x2^2)^(x1^2 + 1);
    startp = [-1.5,1.25]';
    xmin = [0,0]';
    fmin = 0;
    domains = {[-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f6_sym()
    syms x1 x2 x3
    f = (x1^2 + x2^3 - x3^4)^2 + (2*x1*x2*x3)^2 + (2*x1*x2-3*x2*x3+x1*x3)^2;
    startp = [10,10,10]';
    xmin = [0.0048605831,0.0016994507, 0]';
    fmin = 0;
    domains = {[-inf,inf], [-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f7_sym()
    syms x1 x2 x3
    f = x1^2 + (x2 + x2^2)^2 + (-1 + exp(x3))^2;
    startp = [2, 3, -8]';
    xmin = [3.02442e-08, -1, 2.40687e-08]';
    fmin = 1.58063e-15;
    domains = {[-inf,inf], [-inf,inf], [-inf,inf]};
end

function [f,startp,xmin,fmin,domains] = f8_sym()
    syms x1 x2 x3 x4
    f = (x1-1)^2 + (x1-sqrt(x2))^2 + (x2-sqrt(x3))^2 + (x3-sqrt(x4))^2;
    startp = [0.1,0.1,0.1,0.1]';
    xmin = [0.999993, 0.999983, 0.999964, 0.999912]';
    fmin = 1.13719e-10;
    domains = {[-inf,inf], [0,inf], [0,inf], [0,inf]}; # exclude zero take into account grad
end


functions_info = {@f1_sym, 
            @f2_sym, 
            @f3_sym, 
            @f4_sym, 
            @f5_sym, 
            @f6_sym, 
            @f7_sym, 
            @f8_sym};

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
                        'hanger_zhang_bisection_goldstein',
                        'grippo_backtracking_armijo',
                        'grippo_bisection_wolfe_weak',
                        'grippo_wolfe_strong',
                        'grippo_bisection_goldstein'};

table = simulation(functions_info, line_search_methods);

dfile = 'simulation1.csv';
if exist(dfile, 'file') ; delete(dfile); end
cell2csv(dfile,table)
