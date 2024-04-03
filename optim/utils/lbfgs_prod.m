function [pk] = lbfgs_prod(grad_x,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag)
    % BFGS Search Direction
    %
    % This function returns the (L-BFGS) approximate inverse Hessian,
    % multiplied by the negative gradient

    % Set up indexing
    [nVars,m] = size(S);
    if lbfgs_start == 1
        ind = 1:lbfgs_end;
        memory = lbfgs_end-lbfgs_start+1;
    else
        ind = [lbfgs_start:m 1:lbfgs_end];
        memory = m;
    end
    al = zeros(memory,1);
    be = zeros(memory,1);

    pk = -grad_x;
    for j = 1:length(ind)
        i = ind(end-j+1);
        al(i) = (S(:,i)'*pk)/YS(i);
        pk = pk-al(i)*Y(:,i);
    end

    % Multiply by Initial Hessian
    pk = Hdiag*pk;

    for i = ind
        be(i) = (Y(:,i)'*pk)/YS(i);
        pk = pk + S(:,i)*(al(i)-be(i));
    end
end