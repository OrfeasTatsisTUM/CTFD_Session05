if s~=0
    %% Sort Data from less to more nodes
    [n_sorted, I] = sort(n);
    RAM_sorted = vertcat(RAM(I).bytes);
    sRAM_sorted = vertcat(sRAM(I).bytes);
    t_sorted = t(I);
    st_sorted = st(I);

    %% Plots
    figure (2)
    plot(n_sorted, RAM_sorted)
    hold on
    plot(n_sorted, sRAM_sorted)
    xlabel('Grid Size'); ylabel('storage [MB]');
    title("RAM spent (backslash solver)");
    legend('full', 'sparce', 'Location','northwest');
    set(gcf, 'Position',[10,150,620,550]);

    figure (3)
    plot(n_sorted, t_sorted)
    hold on
    plot(n_sorted, st_sorted)
    xlabel('Grid Size'); ylabel('time [sec]');
    title("Time spent (backslash solver)");
    legend('full', 'sparce', 'Location','northwest');
    set(gcf, 'Position',[640,150,620,550]);
else

    if strcmp(solution, 'Test')
        if z == 1 
            %% Random Diagonal Dominant matrix

            figure(i)
            semilogy(itVec,resVec, 'r-');
            title("Random Diagonal Dominant matrix | " + method + " Method")
            xlabel('Iteration')
            ylabel('Residual norm')
            grid on

            % Inability to converge
            if size(itVec,2) == max_iter
                fprintf(2, 'For %s method with Random Diagonal Dominant matrix and tolerance: ', method)
                fprintf(2, "%s the program can't converge\n", tol)
                text(max(xlim), min(ylim),["Didn't converge"], ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'Color','r')
            end
            if ro == 1; text(0.5, 0.5, ["Spectral Radius bigger than 1"], 'Color','r','HorizontalAlignment','center'); end

            % Plot position
            if strcmp(method, 'Jacobi'); set(gcf, 'Position',[10,350,620,550]);
            elseif strcmp(method, 'GaussSeidel'); set(gcf, 'Position',[10,300,620,550]);
            elseif strcmp(method, 'SOR'); set(gcf, 'Position',[10,250,620,550]);
                text(max(xlim), max(ylim),["Relaxation: " + num2str(relax)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
            end
        
        elseif z == 2
            %% Random Three Diagonal Dominant matrix

            figure(i)
            semilogy(itVec,resVec, 'r-');
            title("Random Three Diagonal Dominant matrix | " + method + " Method")
            xlabel('Iteration')
            ylabel('Residual norm')
            grid on

            % Inability to converge
            if size(itVec,2) == max_iter
                fprintf(2, 'For %s method Random Three Diagonal Dominant matrix and tolerance: ', method)
                fprintf(2, "%s the program can't converge\n", tol)
                text(max(xlim), min(ylim),["Didn't converge"], ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'Color','r')
            end
            if ro == 1; text(0.5, 0.5, ["Spectral Radius bigger than 1"], 'Color','r','HorizontalAlignment','center'); end

            % Plot position
            if strcmp(method, 'Jacobi'); set(gcf, 'Position',[640,350,620,550]);
            elseif strcmp(method, 'GaussSeidel'); set(gcf, 'Position',[640,300,620,550]);
            elseif strcmp(method, 'SOR'); set(gcf, 'Position',[640,250,620,550]);
                text(max(xlim), max(ylim),["Relaxation: " + num2str(relax)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
            end

        else
            %% Fully Random matrix

            figure(i)
            semilogy(itVec,resVec, 'r-');
            title("Fully Random matrix | " + method + " Method")
            xlabel('Iteration')
            ylabel('Residual norm')
            grid on

            % Inability to converge
            if size(itVec,2) == max_iter
                fprintf(2, 'For %s method with Fully Random matrix and tolerance: ', method)
                fprintf(2, "%s the program can't converge\n", tol)
                text(max(xlim), min(ylim),["Didn't converge"], ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'Color','r')
            end
            if ro == 1; text(0.5, 0.5, ["Spectral Radius bigger than 1"], 'Color','r','HorizontalAlignment','center'); end
            
            % Plot position
            if strcmp(method, 'Jacobi'); set(gcf, 'Position',[1270,350,620,550]);
            elseif strcmp(method, 'GaussSeidel'); set(gcf, 'Position',[1270,300,620,550]);
            elseif strcmp(method, 'SOR'); set(gcf, 'Position',[1270,250,620,550]);
                text(max(xlim), max(ylim),["Relaxation: " + num2str(relax)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
            end


        end
    else
        figure(i)
        semilogy(itVec,resVec, 'r-');
        title("Discretization of 2D FVM | " + method + " Method")
        xlabel('Iteration')
        ylabel('Residual norm')
        grid on

        % Inability to converge
        if size(itVec,2) == max_iter
            fprintf(2, 'For %s method with discretized matrix (with 2D-FVM) and tolerance: ', method)
            fprintf(2, "%s the program can't converge\n", tol)
            text(max(xlim), min(ylim),["Didn't converge"], ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'Color','r')
        end
        if ro == 1; text(0.5, 0.5, ["Spectral Radius bigger than 1"], 'Color','r','HorizontalAlignment','center'); end

        % Plot position
        if strcmp(method, 'Jacobi'); set(gcf, 'Position',[10,200,620,550]);
        elseif strcmp(method, 'GaussSeidel'); set(gcf, 'Position',[640,200,620,550]);
        elseif strcmp(method, 'SOR'); set(gcf, 'Position',[1270,200,620,550]);
            text(max(xlim), max(ylim),["Relaxation: " + num2str(relax)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
        end
    end
end
