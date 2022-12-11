clear; close all; clc;

%% Orfeas Emmanouil, Tatsis
%% Fernando, Cruz Ceravalls
%% Yuechen, Chen

%% SESSION_05
%  TUM - Ass. Professorship for Thermo Fluid Dynamics
%  WS022-023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady and unsteady heat equation in a 
% non-Cartesian Grid by the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check RAM & solution time (Session 05)
s = 1;
for dimX = 5:11:38
    for dimY = 5:10:35

% Initialize variables
        InitFVM

% Set up the mesh
        [X, Y] = setUpMesh(dimY, dimX, l, formfunction);

% Fill matrix A and vector B. Solve the linear system.
        [T, t(s), st(s), RAM(s), sRAM(s), n(s), A]   = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, TimeIntegrType, theta, simulationType, s);
        if s == 1
            figure(1)
            spy(A,'ro',2)
            set(gcf, 'Position',[1270,150,550,550])
        end

        s = s + 1;
    end
end

%% Use iterative solvers
s = 0;

% Initialize variables
InitFVM

% Set up the mesh
[X, Y] = setUpMesh(dimY, dimX, l, formfunction);

% Fill matrix A and vector B. Solve the linear system.
[T, t(s), st(s), RAM(s), sRAM(s), n(s), A]   = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf, dt, tend, TimeIntegrType, theta, simulationType, s);

%% Make some plots
% postprocess;
postprocess_Session05;