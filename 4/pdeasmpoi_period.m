function u=pdeasmpoi_period(p,e,t)
% Slightlt modified from "pdeasmpoi.m" from the handouts

it1=t(1,:); it2=t(2,:); it3=t(3,:); 
np=size(p,2);
[ar,g1x,g1y,g2x,g2y,g3x,g3y]=pdetrg(p,t);

% Period-Boudary on Segment 2 & 4
p2=[[];[]]; p4=[[];[]];
for j=1:np
    %if(abs(p(2,j))<1e-6 || abs(p(2,j)-1)<1e-6)
    %    continue
    %end
    if(abs(p(1,j)-1)<1e-6)
        p2=[p2 [j;p(2,j)]];
    end
    if(abs(p(1,j))<1e-6)
        p4=[p4 [j;p(2,j)]];
    end
end
p2=(sortrows(p2',2))';
p4=(sortrows(p4',2))';
for j=1:size(p2,2)
    it1(it1==p2(1,j))=p4(1,j);
    it2(it2==p2(1,j))=p4(1,j);
    it3(it3==p2(1,j))=p4(1,j);
end

c3=(g1x.*g2x+g1y.*g2y).*ar;
c1=(g2x.*g3x+g2y.*g3y).*ar;
c2=(g3x.*g1x+g3y.*g1y).*ar;

A=sparse(it1,it2,c3,np,np);
A=A+sparse(it2,it3,c1,np,np);
A=A+sparse(it3,it1,c2,np,np);
A=A+A.';
A=A+sparse(it1,it1,-c2-c3,np,np);
A=A+sparse(it2,it2,-c3-c1,np,np);
A=A+sparse(it3,it3,-c1-c2,np,np);

f=ar/3;
F=sparse(it1,1,f,np,1);
F=F+sparse(it2,1,f,np,1);
F=F+sparse(it3,1,f,np,1);

% Dirichlet-0 on Segment 1 & 3
ie=zeros(1,np);
for j=1:size(e,2)
    if(e(5,j)==1 || e(5,j)==3)
        ie(e(1,j))=1;
        ie(e(2,j))=1;
    end
end
ie(p2(1,:))=1; % Important!!! Otherwise A will be singular.
ie=find(ie);
B=speye(np);
B(:,ie)=[];
A=B'*A*B;
F=B'*F;

un=A\F;
u=B*un; % Redo the Dirichlet-Boundary
u(p2(1,:))=u(p4(1,:)); % Redo the Period-Boudary modification\
u=full(u);
