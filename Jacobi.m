%% Script written by: Fernando, Cruz Ceravalls

function [x1, itVec, resVec, ro, t]=Jacobi(A, B, tol, max_iter)
% arrays for plotting res vs. iteration
resVec = [];
itVec = [];

[~,m]=size(B');
s=size(A,1);

step=0;
x0=zeros(size(B'));
x1=x0;
M=diag(A);
N=A-diag(M);

T = -diag(1./M) * N;
tic
if max(abs(eig(T))) < 1  % Spectral Radius
    ro = 0;
    while norm(B' - A*x1)/norm(B') > tol   &&  step < max_iter
        x1=(B'-N*x0)./M;
        x0=x1;
        step=step+1;
        itVec = [itVec step];
        resVec = [resVec norm(B' - A*x1)/norm(B')];
    end
else
    fprintf(2, 'Spectral radius is bigger than 1!\n')
    ro = 1;
end

t=toc;
    if (step < max_iter) && (max(abs(eig(T))) < 1)
        fprintf('Succesfully calculated! (Solution time: %s [s])\n', num2str(t))
    end
end

