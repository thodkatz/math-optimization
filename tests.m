format longE

# setup the path to include the 'utils' directory
directory = fileparts(which(mfilename));
addpath(genpath(directory));


function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end

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
xmins_newton = {[9.999999999999893e-01;9.999999999999787e-01], 
                [9.999999999816617e-01;9.999999999576692e-01], 
                [9.999999547008681e-01;9.999999073667162e-01],
                [1;1],
                [1;1],
                [1;1],
                [1;1]};

fmins_newton = {1.135959703518257e-28, 
                3.533228612015991e-21, 
                2.466142820222457e-15,
                3.534107889272943e-21,
                2.337282696009845e-26,
                2.337282696009845e-26,
                2.466142972801543e-15};

iters_newton = {6, 
                25, 
                58,
                25,
                8,
                8,
                58};

xmins_steepest = {[1;1], 
                [7.314832480328023e-01;5.367666412124242e-01], 
                [6.826021505276559e-01;4.636937784980137e-01],
                [6.830308556342202e-01;4.645442633745788e-01],
                [1.360838364594753e+00;2.372879081831821e+00],
                [2.422295799273501e+00;6.265632381799888e+00],
                [1.395044815354599e+00;8.002568421122733e+00]};

fmins_steepest = {7.238987188855454e-02, 
                7.238987188855454e-02, 
                1.012485080504488e-01,
                1.008642102262282e-01,
                2.727409876821381e+01,
                1.787251590602936e+01,
                3.668176424944542e+03};

iters_steepest = {100, 
                100, 
                100,
                100,
                100,
                100,
                100};

tol = 1e-6;
for i=1:numel(line_search_methods)
    fprintf("\n")
    line_search_method = line_search_methods{i}
    if strcmp(line_search_method, 'wolfe_strong') || strcmp(line_search_method, 'nonmonotone_wolfe_strong')
        c = [1e-4 0.9];
        rho = 2;
    else
        c = 0.1;
        rho = 0.5;
    end

    tic
    [xmin, fmin, iter] = newton(rosen_sym, 
                        [-1.8,-1.8]', 
                        line_search_method, 
                        search_x, 
                        search_y, 
                        c, 
                        rho, 
                        a=1, 
                        eps=1e-6, 
                        max_iters=100, 
                        to_plot=to_plot);
    assert(xmin, xmins_newton{i}, tol)
    assert(fmin, fmins_newton{i}, tol)
    assert(iter, iters_newton{i})
    toc

    if strcmp(line_search_method, "none")
        continue # steepest descent is gonna fail to converge
    end

    fprintf("\n")
    tic
    [xmin, fmin, iter] = steepest_descent(rosen_sym, 
                        [-1.8,-1.8]', 
                        line_search_method, 
                        search_x, 
                        search_y, 
                        c, 
                        rho, 
                        a=1, 
                        eps=1e-6, 
                        max_iters=100, 
                        to_plot=to_plot);
    assert(xmin, xmins_steepest{i}, tol)
    assert(fmin, fmins_steepest{i}, tol)
    assert(iter, iters_steepest{i})
    toc
end