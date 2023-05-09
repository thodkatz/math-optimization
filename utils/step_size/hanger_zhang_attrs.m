function [ck, qk] = hanger_zhang_attrs(f_value, ck, qk)
    # controls the degree of nonmonotonicity
    # 0-1,0: monotone, 1: nonmonotone
    hetak = 1;
    ck = (hetak*qk*ck + f_value)/(qk*hetak + 1);
    qk = qk*hetak + 1;
end