format longE

# setup the path to include the 'utils' directory
directory = fileparts(which(mfilename)); 
addpath(genpath(directory));    


function f = rosen_sym()
    syms x y
    f = 100*(y-x^2)^2 + (1-x)^2;
end

multidim_grad(rosen_sym(), [1.2,1.2]')
multidim_hessian(rosen_sym(), [1.2,1.2]')
