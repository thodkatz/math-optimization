function plot_line_search1(f, xall, steps_all)
    yall = [];
    iterations = 1:numel(xall);
    for i = iterations
        yall = [yall f(xall{i})];
    end

    # we can't plot more than 2 dimensions
    more_2_dims = !is_vector(xall{1});
    if more_2_dims
        rows = 1;
    else
        rows = 3;
    end
    plot_num = 1;
    if not(more_2_dims)
        subplot(rows, 1, plot_num)
        plot_num += 1;
        # search domain ??
        x = -2:0.5:2;
        y = -2:0.5:2;
        % x = linspace(-2,2);
        % y = linspace(-2,2);
        z = zeros(numel(x),numel(y));
        for row=1:numel(x)
            for col=1:numel(y)
                z(row,col) = f([x(row) y(col)]');
            end
        end
        contour(x,y,z,'Fill','on')
        colorbar
        hold on
        x1all = [];
        x2all = [];
        for i=iterations
            x1 = xall{i}(1);
            x2 = xall{i}(2);
            plot(x1,x2,'h','color','red','MarkerSize',8)
            x1all = [x1all x1];
            x2all = [x2all x2];
        end
        subplot(rows, 1, plot_num)
        plot3(x1all,x2all,yall)
        plot_num += 1;
    end
    subplot(rows, 1, plot_num)
    plot(iterations, steps_all, '-o')
    xlabel("Iterations")
    ylabel("alpha")
end