format longE

# setup the path to include the 'utils' directory
directory = fileparts(which(mfilename));
addpath(genpath(directory));


function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end

tol = 1e-13;
assert(multidim_grad(rosen_sym(), [1.2,1.2]'), [115.6; -48], tol)
assert(multidim_hessian(rosen_sym(), [1.2,1.2]'), [1250, -480; -480, 200])


to_plot=false;
search_x = -1.8:0.1:1.2;
search_y = -1.8:0.1:1.2;

line_search_methods = {'none', 
                        'backtracking_armijo', 
                        'wolfe_strong', 
                        'backtracking_wolfe_weak',
                        'nonmonotone_backtracking_armijo', 
                        'nonmonotone_backtracking_wolfe_weak',
                        'nonmonotone_wolfe_strong'};
xmins = {[1;1], 
        [1;1], 
        [1;1],
        [1;1],
        [1;1],
        [1;1],
        [1;1]};

fmins = {5.816493321323984e-28, 
        3.534107889272943e-21, 
        2.466142972801543e-15,
        3.534107889272943e-21,
        2.337282696009845e-26,
        2.337282696009845e-26,
        2.466142972801543e-15};

iters = {6, 
        25, 
        58,
        25,
        8,
        8,
        58};

count = 1;
tol = 1e-6;
for i=1:numel(line_search_methods)
    tic
    fprintf("\n")
    line_search_method = line_search_methods{i}
    if strcmp(line_search_method, 'wolfe_strong') || strcmp(line_search_method, 'nonmonotone_wolfe_strong')
        c = [1e-4 0.9];
        rho = 2;
    else
        c = 0.1;
        rho = 0.5;
    end
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
    assert(xmin, xmins{count}, tol)
    assert(fmin, fmins{count})
    assert(iter, iters{count})
    count += 1;
    toc
end