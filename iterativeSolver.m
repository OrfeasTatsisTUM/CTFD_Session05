function T = iterativeSolver(solution, s, A, B, tol, relax, max_iter)

i=4; % count of figures
ro = 0;  % ro=0 -> specral radius <= 1,  ro=1 -> specral radius > 1

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Test (Stage 3) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if strcmp(solution, 'Test')
    %% Random A matrix generation
    
    n = 100; % dimension of Atest

    % Random diagonal dominant matrix
    x = 0.3;                                     % fraction of zeros in off-diagonal
    k = round(n*(n-1)*x);                        % number of zeros in off-diagonal

    data = rand(n*(n-1)-k,1);                    % random numbers
    data = [data;zeros(k,1)];                    % the k zeros
    data = data(randperm(length(data)));         % shuffle

    diag_index = 1:n+1:n*n;                      % linear index to all diagonal elements
    offd_index = setdiff(1:n*n,diag_index);      % linear index to all other elements
    Atest1 = zeros(n,n);
    Atest1(offd_index) = data;                   % set off-diagonal elements to data
    Atest1(diag_index) = abs(sum(Atest1,1))+0.1; % set diagonal elements to sum of columns + 0.1
    Atest(:,:,1) = Atest1';                      % columns->rows

    % Random three-diagonal dominant matrix
    Atest2 = diag(rand(1,n-1),1) + diag(rand(1,n-1),-1);  % one diagonal above & one below -> random 
    Atest2(diag_index) = abs(sum(Atest2,1))+0.1;          % set diagonal elements to sum of columns + 0.1
    Atest(:,:,2) = Atest2';                               % columns->rows

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

