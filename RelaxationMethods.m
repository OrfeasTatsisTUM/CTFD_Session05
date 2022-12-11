function Practica1
 clc

A = [-4, 1, 0, 1, 0, 0, 0, 0, 0;
     1, -4, 1, 0, 1, 0, 0, 0, 0;
     0, 1, -4, 0, 0, 1, 0, 0, 0;
     1, 0, 0, -4, 1, 0, 1, 0, 0;
     0, 1, 0, 1, -4, 1, 0, 1, 0;
     0, 0, 1, 0, 1, -4, 0, 0, 1;
     0, 0, 0, 1, 0, 0, -4, 1, 0;
     0, 0, 0, 0, 1, 0, 1, -4, 1;
     0, 0, 0, 0, 0, 1, 0, 1, -4];

          
y = [-120; 0; -30; -70; 0; -20; -290; -170; -160];
tolerancia = 0.01;
w=1.2;



[x0,x1,i] = Jacobi(A,y,tolerancia)


end

function r=IsPredomDiag(f)
    [n m]=size(f);
    if m~=n
        r=0;
    else
        r=all(2*diag(abs(f))-sum(abs(f),2)>=-1e-15);
    end      

end

function [x0,x1,i]=Jacobi(A,y,tolerancia)
    [n,m]=size(y);
    s=size(A,1);
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
    while (e>tolerancia)
            x1=(y-N*x0)./M;
            e=norm(x1-x0);
            x0=x1;
            i=i+1;   
    end
end

function [x0,x1,step]=GaussSeidel(A,y,tolerancia)
 [n,m]=size(y);
 step = 0;
 s=size(A,1);
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
    while (norm(x1-x0)>tolerancia)

        x0=x1;
        for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))/A(i,i);
        end
        i=i+1;
     step = step+1;
    end
end

function [x0,x1,step]=Relajacion(A,y,tolerancia,w)
 [n,m]=size(y);
 step = 0;
 s=size(A,1);
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
    while (norm(x1-x0)>tolerancia)
        x0=x1;
        for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i)+(1-w)*x0(i,:);
        end
        i=i+1;
        step = step+1;
    end
end










