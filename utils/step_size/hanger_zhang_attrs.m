function [ck, qk] = hanger_zhang_attrs(f_value, ck, qk)
    % Zhang, H. and Hager, W.W. (2004) A Nonmonotone Line Search Technique and Its
    % Application to Unconstrained Optimization. SIAM Journal on Optimization, 14,
    % 1043- 1056
    %
    % controls the degree of nonmonotonicity
    % 0-1,0: monotone, 1: nonmonotone
    hetak = 1;
    ck = (hetak*qk*ck + f_value)/(qk*hetak + 1);
    qk = qk*hetak + 1;
end