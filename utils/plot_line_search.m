function plot_line_search(f, xall, steps_all)
    yall = [];
    iterations = 1:numel(xall);
    for i = iterations
        yall = [yall f(xall{i})];
    end

    # we can't plot more than 2 dimensions
    more_2_dims = !is_vector(xall{1});
    if more_2_dims
        rows = 2;
    else
        rows = 3;
    end
    plot_num = 1;
    if !more_2_dims
        subplot(rows, 1, plot_num)
        plot_num += 1;
        x1 = [];
        x2 = [];
        for i=iterations
            x1 = [x1 xall{i}(1)];
            x2 = [x2 xall{i}(2)];
        end
        plot(iterations, x1, '-o')
        hold on
        plot(iterations, x2, '-o')
        legend("x1", "x2")
        xlabel("Iterations")
        ylabel("Coordinates")
    end
    subplot(rows, 1, plot_num)
    plot(iterations, yall, '-o')
    xlabel("Iterations")
    ylabel("f(x)")
    subplot(rows, 1, plot_num + 1)
    plot(iterations, steps_all, '-o')
    xlabel("Iterations")
    ylabel("alpha")
end