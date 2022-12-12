function [x1, itVec, resVec, ro]=SOR(A, B, tol, relax, max_iter)
% arrays for plotting res vs. iteration
resVec = [];
itVec = [];

[n,~]=size(B');
step = 0;
s=size(A,1);

x0=zeros(size(B'));
x1=x0;

if max(eig(diag(A) + tril(A, -1))) < 1  % Spectral Radius
    ro = 0;
    for i=1:n
        x1(i,:)=(B(i)'-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*relax/A(i,i);
    end

    while (norm(B' - A*x1)/norm(B') > tol  &&  step < max_iter)
        x0=x1;
        for i=1:n
            x1(i,:)=(B(i)'-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*relax/A(i,i)+(1-relax)*x0(i,:);
        end
        step = step+1;
        itVec = [itVec step];
        resVec = [resVec norm(B' - A*x1)];
    end
else
    fprintf(2, 'Spectral radius is bigger than 1\n')
    ro = 1;
end
end