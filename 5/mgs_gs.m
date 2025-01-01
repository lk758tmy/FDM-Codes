function x=mgs_gs(A,r,x0,m,flag)
% (Local) Gauss-Seidel smoother for multigrid method.
% Copied from the textbook.
x=x0;
if flag==1
    L=tril(A);
    U=triu(A,1);
    for k=1:m
        x=L\(r-U*x);
    end
end
if flag==-1
    L=tril(A,-1);
    U=triu(A);
    for k=1:m
        x=U\(r-L*x);
    end
end