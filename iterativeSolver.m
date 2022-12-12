function T = iterativeSolver(solution, s, A, B, tol, relax, max_iter)

i=4; % count of figures
ro = 0;  % ro=0 -> specral radius <= 1,  ro=1 -> specral radius > 1

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Test (Stage 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strcmp(solution, 'Test')
    %% Random A matrix generation
    n = 100; % dimension of Atest

    % Random diagonal dominant matrix
    Atest = rand(n);
    Atest(logical(eye(size(Atest)))) = n*n;
    Atest(:,:,1) = Atest;      % A should be diagonal dominant GS and J converge, SOR convergence not guaranteed

    % Random three-diagonal dominant matrix
    A1 = diag(rand(1,n-1),1) + diag(rand(1,n-1),-1) + diag(n*n*ones(1,n));
    Atest(:,:,2) = A1;         % A diagonal dominant => GS and J converge, SOR convergence not guaranteed

    % Fully random matrix
    Atest(:,:,3) = rand(n);

    Btest = rand(1,n);

    %% Solving
    
    for z=1:3

        %% Check if conditions are correct
        [n,~]=size(Btest');
        sizeA=size(Atest(:,:,z), 1);
        if ~IsPredomDiag(Atest(:,:,z))
            if z==1
                fprintf(2,"Random Diagonal Dominant matrix is not ")
            elseif z==2
                fprintf(2,"Random Three Diagonal Dominant matrix is not ")
            else
                fprintf(2,"Fully random Matrix is not ")
            end
            fprintf(2,"diagonaly predominant\n")
            break;
        end

        if sizeA~=n
            error('Dimension Error');
        end

        % check if A is singular
        if det(Atest(:,:,z)) == 0
            error('A is singular! :( \n')
        end

        %% Jacobi

        method = 'Jacobi';
        [T, itVec, resVec, ro] = Jacobi(Atest(:,:,z), Btest, tol, max_iter);
        postprocess_Session05;
        i = i+1;

        %% Gauss-Seidel

        method = 'GaussSeidel';
        [T, itVec, resVec, ro] = SOR(Atest(:,:,z), Btest, tol, 1, max_iter);  % Gaus-Seidel method is SOR with relaxation = 1
        postprocess_Session05;
        i = i+1;

        %% Succesive Over Relaxation (SOR)

        method = 'SOR';
        [T, itVec, resVec, ro] = SOR(Atest(:,:,z), Btest, tol, relax, max_iter);
        postprocess_Session05;
        i = i+1;

    end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2D FVM (Stage 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
elseif strcmp(solution, '2D_FVM')
    tol = 1.0e-2;  % Tolerance for solver = '2D_FVM'  (predefined in Scriptum)

    method = 'Jacobi';
    [T, itVec, resVec, ro] = Jacobi(A, B, tol, max_iter);  % Jacobi
    postprocess_Session05;
    i = i+1;
    
    method = 'GaussSeidel';
    [T, itVec, resVec, ro] = SOR(A, B, tol, 1, max_iter);  % Gaus-Seidel (SOR with relaxation = 1)
    postprocess_Session05;
    i = i+1;

    method = 'SOR';
    [T, itVec, resVec, ro] = SOR(A, B, tol, relax, max_iter);  % Succesive Over Relaxation (SOR)
    postprocess_Session05;
    i = i+1;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Incorrect Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else
    error('Incorrect input in variable: solution')
end
end

%% Diagonal Predominant Check Function
function r=IsPredomDiag(f)
[n, m]=size(f);
if m~=n
    r=0;
else
    r=all(2*diag(abs(f))-sum(abs(f),2)>=-1e-15);
end
end

