clear; close all; clc;

%% Orfeas Emmanouil, Tatsis
%% Fernando, Cruz Ceravalls
%% Yuechen, Chen

%% SESSION_05
%  TUM - Ass. Professorship for Thermo Fluid Dynamics
%  WS022-023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid 
% by the Finite Volumes Method using self made solvers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Question = input("Which Stage would you like to solve?\n Press 1: Stage 2\n Press 2: Stage 3, 4 & 5\n Press 3: Stage 6\n");
FLAG = false;
while ~FLAG
if Question == 1
    Stg2 = 1; Stg6 = 0; FLAG = true;
elseif Question == 2
    Stg2 = 0; Stg6 = 0; FLAG = true;
elseif Question == 3
    Stg2 = 0; Stg6 = 1; FLAG = true;
else
    clc; fprintf(2,"Wrong input! Please chose between 1, 2 & 3\n")
    Question = input("Which Stage would you like to solve?\n Press 1: Stage 2\n Press 2: Stage 3, 4 & 5\n Press 3: Stage 6\n");
end
end
%% Check RAM & solution time  with backslash solver for various grid sizes (Session 05 - Stage 2)

if Stg2 ~= 0
    for dimX = 5:11:38
        for dimY = 5:10:35

            % Initialize variables
            InitFVM

            % Set up the mesh
            [X, Y] = setUpMesh(dimY, dimX, l, formfunction);

            % Fill matrix A and vector B. Solve the linear system.
            [T, t(Stg2), st(Stg2), RAM(Stg2), sRAM(Stg2), n(Stg2), A]   = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, TimeIntegrType, theta, simulationType, Stg2, tol, max_iter, relax, solution, Stg6);
            if Stg2 == 1
                figure(1)
                spy(A,'ro',2)
                set(gcf, 'Position',[1270,150,550,550])
                title("Matrix discretized with 2D FVM from grid: " + dimX + "x" + dimY + " ")
            end

            Stg2 = Stg2 + 1;
        end
    end

    % Make some plots
    postprocess_Session05;

else
    %% Use iterative solvers (Session 05 - Stage 3, 4, 5)

    if Stg6 == 0
        % Initialize variables
        InitFVM

        % Set up the mesh
        [X, Y] = setUpMesh(dimY, dimX, l, formfunction);

        % Fill matrix A and vector B. Solve the linear system.
        T   = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, TimeIntegrType, theta, simulationType, Stg2, tol, max_iter, relax, solution, Stg6);

        % Make some plots
        if strcmp(solution, '2D_FVM');  postprocess;  end

    else
    %% Best iterative solver vs gmres (Session 05 - Stage 6)
        for Stg6 = 1:3

            dimX = 21 * Stg6;
            dimY = dimX;

            if Stg6 == 1
                fprintf("Medium size matrix (Node number: %s)\n", num2str(dimX*dimY))
            elseif Stg6 == 2
                fprintf("Large size matrix (Node number: %s)\n", num2str(dimX*dimY))
            else
                fprintf("Very large size matrix (Node number: %s)\n", num2str(dimX*dimY))
            end

            % Initialize variables
            InitFVM

            % Set up the mesh
            [X, Y] = setUpMesh(dimY, dimX, l, formfunction);

            % Fill matrix A and vector B. Solve the linear system.
            T   = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, TimeIntegrType, theta, simulationType, Stg2, tol, max_iter, relax, solution, Stg6);
        end

    end
end