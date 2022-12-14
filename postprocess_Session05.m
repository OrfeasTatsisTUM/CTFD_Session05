%% Script written by: Orfeas-Emmanouil Tatsis

if Stg2 ~= 0 && Stg6 == 0  % Stage 2
    %% Sort Data from less to more nodes
    [n_sorted, I] = sort(n);
    RAM_sorted = vertcat(RAM(I).bytes);
    sRAM_sorted = vertcat(sRAM(I).bytes);
    t_sorted = t(I);
    st_sorted = st(I);

    %% Plots (Stage 2)
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
    
elseif Stg2 == 0 && Stg6 == 0  % Stage 3, 4, 5
    %% Stage 3
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
                fprintf(2, "The program couldn't converge with tolerance %s\n", num2str(tol))
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
                fprintf(2, "The program couldn't converge with tolerance %s\n", num2str(tol))
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
                fprintf(2, "The program couldn't converge with tolerance %s\n", num2str(tol))
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
        %% FVM matrix (Stage 4, 5)
        figure(i)
        semilogy(itVec,resVec, 'r-');
        title("Discretization of 2D FVM | " + method + " Method")
        xlabel('Iteration')
        ylabel('Residual norm')
        grid on

        % Inability to converge
        if size(itVec,2) == max_iter
            fprintf(2, "The program couldn't converge with tolerance %s\n", num2str(tol))
            text(max(xlim), min(ylim),["Didn't converge"], ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'Color','r')
        end
        if ro == 1; text(0.5, 0.5, ["Spectral Radius bigger than 1"], 'Color','r','HorizontalAlignment','center'); end

        % Plot position
        if strcmp(method, 'Jacobi'); set(gcf, 'Position',[10,200,620,550]);
            text(max(xlim), max(ylim),["Solving time: " + num2str(tJac) + " [s]"], 'VerticalAlignment','top', 'HorizontalAlignment','right')
        elseif strcmp(method, 'GaussSeidel'); set(gcf, 'Position',[640,200,620,550]);
            text(max(xlim), max(ylim),["Solving time: " + num2str(tGS) + " [s]"], 'VerticalAlignment','top', 'HorizontalAlignment','right')
        elseif strcmp(method, 'SOR'); set(gcf, 'Position',[1270,200,620,550]);
            text(max(xlim), max(ylim),["Relaxation: " + num2str(relax)], 'VerticalAlignment','top', 'HorizontalAlignment','right')
            text(max(xlim), 0.5*max(ylim),["Solving time: " + num2str(tSOR) + " [s]"], 'VerticalAlignment','top', 'HorizontalAlignment','right')
        end
    end

elseif Stg6 ~=0  % Stage 6

    if ~strcmp(method, 'gmres')
    %% Best self made iterative solver
        i = Stg6+19;
        figure(i)
        semilogy(itVec,resVec, 'r-');
        if Stg6 == 1; title("Medium size 2D-FVM | " + method + " Method");
        elseif Stg6 == 2; title("Large size 2D-FVM | " + method + " Method");
        else; title("Very large size 2D-FVM | " + method + " Method");
        end

        xlabel('Iteration')
        ylabel('Residual norm')
        grid on
        set(gcf, 'Position',[300,200,620,550]);
        text(max(xlim), 0.5*max(ylim),["Solving time: " + num2str(t_own) + " [s]"], 'VerticalAlignment','top', 'HorizontalAlignment','right')
        if strcmp(method, 'SOR'); text(max(xlim), max(ylim),["Relaxation: " + num2str(relax)], 'VerticalAlignment','top', 'HorizontalAlignment','right'); end

        % Inability to converge
        if size(itVec,2) == max_iter
            fprintf(2, "Iterative solver couldn't converge with tolerance %s\n", num2str(tol))
            text(max(xlim), min(ylim),["Didn't converge"], ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','right', 'Color','r')
        end
        if ro == 1; text(0.5, 0.5, ["Spectral Radius bigger than 1"], 'Color','r','HorizontalAlignment','center'); end

    else
    %% gmres
        i = Stg6+22;
        figure(i)
        semilogy(0:length(resVec)-1,resVec/norm(B), 'r-');
        if Stg6 == 1; title("Medium size 2D-FVM | " + method + " Method");
        elseif Stg6 == 2; title("Large size 2D-FVM | " + method + " Method");
        else; title("Very large size 2D-FVM | " + method + " Method");
        end
        xlabel('Iteration')
        ylabel('Residual norm')
        grid on
        set(gcf, 'Position',[980,200,620,550]);
        text(max(xlim), 0.5*max(ylim),["Solving time: " + num2str(t) + " [s]"], 'VerticalAlignment','top', 'HorizontalAlignment','right')
        if conv
            fprintf(2, "Didn't converge with tolerance %s\n", num2str(tol))
        else
            fprintf('Succesfully calculated! (Solution time: %s [s])\n', num2str(t))
        end
    end
end
