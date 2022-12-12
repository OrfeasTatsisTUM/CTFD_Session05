function [x1, itVec, resVec]=SOR(A, y, tol, w, maxiter)
% arrays for plotting res vs. iteration
resVec = [];
itVec = [];

[n,~]=size(y);
step = 0;
s=size(A,1);

x0=zeros(size(y));
x1=x0;

for i=1:n
    x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i);
end

while (norm(x1-x0)>tol && step < maxiter)
    x0=x1;
    for i=1:n
        x1(i,:)=(y(i,:)-A(i,1:i-1)*x1(1:i-1,:)-A(i,i+1:n)*x0(i+1:n,:))*w/A(i,i)+(1-w)*x0(i,:);
    end
    step = step+1;
    itVec = [itVec step];
    resVec = [resVec norm(y - A*x1)];
end
end