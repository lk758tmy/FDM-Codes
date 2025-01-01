function [u,steps]=mg_vcycle(A,F,I,u0,m,k,tol)
% Multigrid V-cycle iteration.
% Copied from the textbook.
if k==1
    u=A\F;
    steps=1;
else
    u=u0;
    r=F-A*u;
    error0=max(abs(r));
    error=error0;
    steps=0;
    while error>tol*error0
        Br=mgp_vcycle(A,r,I,m,k); % Multigrid precondtioner
        u=u+Br;
        r=F-A*u;
        error=max(abs(r));
        steps=steps+1;
    end
end
