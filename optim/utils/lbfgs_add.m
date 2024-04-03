function [S,Y,YS,lbfgs_start,lbfgs_end,Hdiag,skipped] = lbfgs_add(y,s,S,Y,YS,lbfgs_start,lbfgs_end,Hdiag)
    ys = y'*s;
    skipped = 0;
    m = size(S,2);
    if ys > 1e-10
        if lbfgs_end < m
            lbfgs_end = lbfgs_end+1;
            if lbfgs_start ~= 1
                if lbfgs_start == m
                    lbfgs_start = 1;
                else
                    lbfgs_start = lbfgs_start+1;
                end
            end
        else
            lbfgs_start = min(2,m);
            lbfgs_end = 1;
        end
        
        S(:,lbfgs_end) = s;
        Y(:,lbfgs_end) = y;
        YS(lbfgs_end) = ys;
        
        % Update scale of initial Hessian approximation
        Hdiag = ys/(y'*y);
    else
        skipped = 1;
    end
end