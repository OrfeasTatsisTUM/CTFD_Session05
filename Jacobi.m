function [x1, itVec, resVec]=Jacobi(A, y, tol, maxiter)
% arrays for plotting res vs. iteration
resVec = [];
itVec = [];

[~,m]=size(y);
s=size(A,1);

step=0;
x0=zeros(size(y));
e=realmax;
M=diag(A);
N=A-diag(M);
M=M*ones(1,m);

while (e>tol  && step < maxiter)
    x1=(y-N*x0)./M;
    e=norm(x1-x0);
    x0=x1;
    step=step+1;
    itVec = [itVec step];
    resVec = [resVec norm(y - A*x1)];
end
end

