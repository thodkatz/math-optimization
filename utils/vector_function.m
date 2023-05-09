function ret = vector_function(fun_handler, vector)
    if nargin(fun_handler) == numel(vector)
        ret = fun_handler(num2cell(vector){:});
    elseif nargin(fun_handler) ~= 0
        arg_list = get_arg_names(fun_handler);
        filt_arg = [];
        for i=1:numel(arg_list)
            arg_ith = arg_list{i};
            filt_arg(end+1) = vector(get_index_of_arg(arg_ith));
        end
        filt_arg = filt_arg';
        ret = fun_handler(num2cell(filt_arg){:});
    else
        ret = fun_handler();
    end
end

function arg_list = get_arg_names(fun_handler)
    name = functions(fun_handler).function;
    arg_pattern = '^\s*@?\(?([^)]+)\)?';
    arg_str = regexp(name,arg_pattern,'tokens','once'){1};
    arg_str = strrep(arg_str,' ','');
    arg_list = strsplit(arg_str,',');
end

function idx = get_index_of_arg(arg)
    % expected format 'x<number>'
    idx = str2num(arg(2:end));
end