function Practica1
 clc

% dimension of A
n = 100;
tol = 1.0e-100;
w = 1.9;
y = rand(n,1);

% arrays for plotting res vs. iteration
resVec = [];
itVec = [];

% create random diagonal dominant matrix 
fprintf('random diagonal dominant matrix \n');

A = rand(n);
% A should be diagonal dominant GS and J converge, SOR convergence not guaranteed
A(logical(eye(size(A)))) = n*n

% create random three-diagonal dominant matrix
fprintf('random three-diagonal dominant matrix \n');
% A diagonal dominant => GS and J converge, SOR convergence not guaranteed
%A = diag(rand(1,n-1),1) + diag(rand(1,n-1),-1) + diag(n*n*ones(1,n))



%% Methods
method='JAC';



%%Jacobi
switch method
    case 'JAC'
    [n,m]=size(y);
    s=size(A,1);
    step = 0;
    
    %Check if conditions are correct
    if ~IsPredomDiag(A)
        error('Not diagonaly predominant');
    end
    if s~=n
        error('Dimmension Error');
    end
    i=0;
    x0=zeros(size(y));
    e=realmax;
    M=diag(A);
    N=A-diag(M);
    M=M*ones(1,m);
    while (e>tol)
            x1=(y-N*x0)./M;
            e=norm(x1-x0);
            x0=x1;
            i=i+1;
            step = step+1;
            
            itVec = [itVec step];
            resVec = [resVec norm(y - A*x1)];
    end

%%GaussSeidel   
    case 'GAU'
     [n,m]=size(y);
     step = 0;
     s=size(A,1);
 
     %Check if conditions are correct
    if ~IsPredomDiag(A)
        error('Not diagonaly predominant');
    end
    if s~=n
        error('Dimmension Error');
    end
    i=0;
    x0=zeros(size(y));
    x1=x0;
    for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))/A(i,i);
    end
    while (norm(x1-x0)>tol)

        x0=x1;
        for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))/A(i,i);
        end
        i=i+1;
        step = step+1;
    end
    
    case SOR
        
     [n,m]=size(y);
     step = 0;
     s=size(A,1);
     %Check if conditions are correct
    if ~IsPredomDiag(A)
        error('Not diagonaly predominant');
    end
    if s~=n
        error('Dimmension Error');
    end
    
    i=0;
    x0=zeros(size(y));
    x1=x0;
    for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i);
    end
    while (norm(x1-x0)>tol)
        x0=x1;
        for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i)+(1-w)*x0(i,:);
        end
        i=i+1;
        step = step+1;
    end
end

step
x1
itVec
resVec


%% plot residual vs. iteration

plot(itVec,resVec, 'r-');
%title('Random diag dominant Matrix, Jacobi Method')
%title('Random diag dominant Matrix, GS Method')
%title('Random diag dominant Matrix, SOR Method, \omega = 1.9')
%title('Random three-diag dominant Matrix, Jacobi Method')
%title('Random three-diag dominant Matrix, GS Method')
%title('Random three-diag dominant Matrix, SOR Method \omega = 1.9')
title('Random full Matrix, Jacobi Method')
xlabel('iteration')
ylabel('residual norm')
grid on
end

function r=IsPredomDiag(f)
    [n m]=size(f);
    if m~=n
        r=0;
    else
        r=all(2*diag(abs(f))-sum(abs(f),2)>=-1e-15);
    end      

end


        
function [x0,x1,i]=Jacobi(A,y,tol)
    [n,m]=size(y);
    s=size(A,1);
    
 %Check if conditions are correct
    if ~IsPredomDiag(A)
        error('Not diagonaly predominant');
    end
    if s~=n
        error('Dimmension Error');
    end
    i=0;
    x0=zeros(size(y));
    e=realmax;
    M=diag(A);
    N=A-diag(M);
    M=M*ones(1,m);
    while (e>tol)
            x1=(y-N*x0)./M;
            e=norm(x1-x0);
            x0=x1;
            i=i+1;   
    end
end
function [x0,x1,step]=GaussSeidel(A,y,tol)
 [n,m]=size(y);
 step = 0;
 s=size(A,1);
 
 %Check if conditions are correct
    if ~IsPredomDiag(A)
        error('Not diagonaly predominant');
    end
    if s~=n
        error('Dimmension Error');
    end
    i=0;
    x0=zeros(size(y));
    x1=x0;
    for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))/A(i,i);
    end
    while (norm(x1-x0)>tol)

        x0=x1;
        for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))/A(i,i);
        end
        i=i+1;
     step = step+1;
    end
end
function [x0,x1,step]=SOR(A,y,tol,w)
 [n,m]=size(y);
 step = 0;
 s=size(A,1);
 %Check if conditions are correct
    if ~IsPredomDiag(A)
        error('Not diagonaly predominant');
    end
    if s~=n
        error('Dimmension Error');
    end
    
    i=0;
    x0=zeros(size(y));
    x1=x0;
    for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i);
    end
    while (norm(x1-x0)>tol)
        x0=x1;
        for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i)+(1-w)*x0(i,:);
        end
        i=i+1;
        step = step+1;
    end
end









