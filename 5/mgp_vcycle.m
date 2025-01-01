function Br=mgp_vcycle(A,r,I,m,k)
% Multigrid V-cycle precondtioner.
% Copied from the textbook.

if k==1
    Br=A\r;
else
    Ik=I{k-1}; % Prolongation matrix
    y=zeros(size(r));
    y=mgs_gs(A,r,y,m,1); % Pre-smoothing
    r1=r-A*y;
    r1=Ik'*r1;
    B=Ik'*A*Ik;
    Br1=mgp_vcycle(B,r1,I,m,k-1);
    y=y+Ik*Br1;
    y=mgs_gs(A,r,y,m,-1); % Post-smoothing
    Br=y;
end