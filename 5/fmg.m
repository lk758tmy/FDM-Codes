function u=fmg(A,F,B,I,l,m,depth)
% 
% "nmg": Number of multigrid iterations.
% "nsm": Number of smoothing iterations.
%
% 根据教材上的代码修改。

u=A{1}\F{1};
for i=2:depth
    u=I{i-1}*u; % Initial value
    for l=1:l % Multigrid iteration
        r=F{i}-A{i}*u;
        Br=mgp_vcycle(A{i},r,I,m,i); % Multigrid precondtioner
        u=u+Br;
    end
end