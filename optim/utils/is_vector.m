function ret = is_vector(point)
    ret = numel(size(point)) == 2 && size(point)(1) > 1 && size(point)(2) == 1;
end