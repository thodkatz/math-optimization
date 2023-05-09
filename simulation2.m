%% Simulation for the functions listed on Chapter 4 of the doc

dfile = 'simulation2.txt';
if exist(dfile, 'file') ; delete(dfile); end
diary(dfile)


format longE

# setup the path to include the 'utils' directory
directory = pwd
addpath(genpath(directory))

function [f,startp,xmin,fmin] = rosen_extended_sym(n)
    syms x1 x2 real
    f = 100*(x2-x1^2)^2 + (1-x1)^2;

    startp = [-1.2,1];
    xmin = [1,1];

    for i = 3:n
        syms(['x', num2str(i)], 'real');
        ith = num2str(i);
        ith_prev = num2str(i-1);
        f = f + (1 - sym(['x', ith_prev]))^2 + 100*(sym(['x', ith]) - sym(['x', ith_prev])^2)^2;

        if mod(i,2) == 0
            startp(end+1) = 1;
        else
            startp(end+1) = -1.2;
        end
        xmin(end+1) = 1;
    end
    startp = startp';
    xmin = xmin';
    fmin = 0;
end

function [f,startp,xmin,fmin] = wood()
    syms x1 x2 x3 x4 real
    f = 100*(x1^2-x2)^2 + (x1-1)^2 + (x3-1)^2 + 90*(x3^2 - x4)^2 + 10.1*((x2-1)^2 + (x4-1)^2) + 19.8*(x2-1)*(x4-1);
    startp = [-3,-1,-3,-1]';
    xmin = [1,1,1,1]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = powell()
    syms x1 x2 x3 x4 real
    f = (x1+10*x2) + 5*(x3-x4)^2 + (x2-2*x3)^4 + 10*(x1-x4)^2;
    startp = [3,-1,0,1]';
    xmin = [0,0,0,0]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = cube()
    syms x1 x2 real
    f = 100*(x2-x1^3)^2 + (1-x1)^2
    startp = [-1.2,1]';
    xmin = [1,1]';
    fmin = 0;
end

function [f,startp,xmin,fmin] = trig(n)
    syms x1
    f = n + (1-cos(x1)) - sin(x1) - cos(x1);
    startp = [1/(5*n)];
    xmin = [0];
    for j=2:n
        syms(['x', num2str(j)], 'real');
        sym_var = sym(['x', num2str(j)]);
        f += cos(sym_var);
        startp(end+1) = 1/(5*n);
        xmin(end+1) = 0;
    end

    for i=2:n
        sym_var = sym(['x', num2str(i)]);
        f += n + i*(1-cos(sym_var)) - sin(sym_var) - cos(x1);
        for j=2:n
            sym_var = sym(['x', num2str(j)]);
            f += cos(sym_var);
        end
    end
    startp = startp';
    xmin = xmin';
    fmin = 0;
end


function [f1,f2,startp,xmin,fmin] = helical()
    syms x1 x2 x3
    theta = atan(x2/x1)/(2*pi);
    f1 = 100*((x3-10*theta)^2 + (sqrt(x1^2+x2^2)-1)^2) + x3^2;

    theta += 1/2;
    f2 = 100*((x3-10*theta)^2 + (sqrt(x1^2+x2^2)-1)^2) + x3^2;
    startp = [-1,0,0]';
    xmin = [1,0,0]';
    fmin = 0;
end

functions = {rosen_extended_sym, 
            wood, 
            powell, 
            cube, 
            trig, 
            helical};

starting_points = {[-1.8, -1.8]',
                    [-3,-3]',
                    [-9.5,9.5]',
                    [1.9,-1.9]',
                    [-1.5,1.25]',
                    [10,10,10]'};
expected_fmin = {0,
                0,
                0,
                -6.4931,
                0,
                0};
expected_xmin = {[1,1]',
                    [0,0]',
                    [1,3]',
                    [-1.1515,0.5455]',
                    [0,0]',
                    [0.0048605831,0.0016994507, 0]'};

to_plot=false;
% redundant search_x and search_y because we don't plot
% todo use nargin, to configure the inputs of the functions better
search_x = -1.8:0.1:1.2;
search_y = -1.8:0.1:1.2;

line_search_methods = {'hanger_zhang_backtracking_armijo'};

diary off