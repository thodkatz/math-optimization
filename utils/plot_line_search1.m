function plot_line_search1(f, yall, xall, steps_all, search_domain_x, search_domain_y)
    iterations = 1:numel(xall);
    # we can't plot more than 2 dimensions
    more_2_dims = numel(xall{1}) > 2;
    if more_2_dims
        rows = 2;
    else
        rows = 6;
    end
    plot_num = 1;
    if not(more_2_dims)
        subplot(rows, 1, plot_num)
        plot_num += 1;
        x = search_domain_x;
        y = search_domain_y;
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
            plot(x1,x2,'o','color','red','MarkerSize',8)
            x1all = [x1all x1];
            x2all = [x2all x2];
        end
        subplot(rows,1,plot_num)
        plot3(x1all,x2all,yall)
        plot_num += 1;

        subplot(rows,1,plot_num)
        plot(iterations,x1all)
        ylabel("x1")
        plot_num += 1;

        subplot(rows,1,plot_num)
        plot(iterations,x2all)
        ylabel("x2")
        plot_num += 1;
    end

    subplot(rows, 1, plot_num)
    plot(iterations,yall)
    ylabel("f")
    plot_num += 1;

    subplot(rows, 1, plot_num)
    plot(iterations, steps_all, '-o')
    xlabel("Iterations")
    ylabel("alpha")
end