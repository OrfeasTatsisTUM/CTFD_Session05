%% Script written by: Orfeas-Emmanouil Tatsis

function T = iterativeSolver(solution, Stg2, A, B, tol, relax, max_iter, Stg6)

if Stg6 == 0
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

            if z==1
                fprintf("Random Diagonal Dominant matrix:\n")
            elseif z==2
                fprintf("Random Three Diagonal Dominant matrix:\n")
            else
                fprintf("Fully random Matrix:\n  ")
            end

            %% Check if conditions are correct
            [n,~]=size(Btest');
            sizeA=size(Atest(:,:,z), 1);
            if ~IsPredomDiag(Atest(:,:,z))
                fprintf(2,"Not diagonaly predominant!\n")
                break;
            end

            % Check dimensions
            if sizeA~=n
                fprintf(2,'Dimension Error!\n');
                break;
            end

            % check if A is singular
            if det(Atest(:,:,z)) == 0
                fprintf(2,'Matrix is singular!\n')
                break;
            end

            %% Jacobi

            fprintf("  Jacobi Method:\t\t\t\t")
            method = 'Jacobi';
            [T, itVec, resVec, ro, tJac] = Jacobi(Atest(:,:,z), Btest, tol, max_iter);
            postprocess_Session05;
            i = i+1;

            %% Gauss-Seidel

            fprintf("  Gauss-Seidel Method:\t\t\t")
            method = 'GaussSeidel';
            [T, itVec, resVec, ro, tGS] = SOR(Atest(:,:,z), Btest, tol, 1, max_iter);  % Gaus-Seidel method is SOR with relaxation = 1
            postprocess_Session05;
            i = i+1;

            %% Succesive Over Relaxation (SOR)

            fprintf("  SOR Method (ω=%s):\t\t\t", num2str(relax))
            method = 'SOR';
            [T, itVec, resVec, ro, tSOR] = SOR(Atest(:,:,z), Btest, tol, relax, max_iter);
            postprocess_Session05;
            i = i+1;

        end

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 2D FVM (Stage 4) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elseif strcmp(solution, '2D_FVM')
        tol = 1.0e-2;  % Tolerance for solver = '2D_FVM'  (predefined in Scriptum)

        fprintf("  Jacobi Method:\t\t\t\t")
        method = 'Jacobi';
        [T, itVec, resVec, ro, tJac] = Jacobi(A, B, tol, max_iter);  % Jacobi
        postprocess_Session05;
        i = i+1;

        fprintf("  Gauss-Seidel Method:\t\t\t")
        method = 'GaussSeidel';
        [T, itVec, resVec, ro, tGS] = SOR(A, B, tol, 1, max_iter);  % Gaus-Seidel (SOR with relaxation = 1)
        postprocess_Session05;
        i = i+1;

        fprintf("  SOR Method (ω=%s):\t\t\t", num2str(relax))
        method = 'SOR';
        [T, itVec, resVec, ro, tSOR] = SOR(A, B, tol, relax, max_iter);  % Succesive Over Relaxation (SOR)
        postprocess_Session05;
        i = i+1;

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Incorrect Input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else
        error('Incorrect input in variable: solution')
    end

else

    fprintf("    SOR Method (ω=%s):\t\t\t", num2str(relax))
    method = 'SOR';
    [T, itVec, resVec, ro, t_own] = SOR(A, B, tol, relax, max_iter);  % Succesive Over Relaxation (SOR)
    postprocess_Session05;

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

