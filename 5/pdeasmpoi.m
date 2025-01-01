function [A,F,B]=pdeasmpoi(p,e,t)

% 初步工作
it1=t(1,:);
it2=t(2,:);
it3=t(3,:);
np=size(p,2);
[ar,g1x,g1y,g2x,g2y,g3x,g3y]=pdetrg(p,t);

% 装配单元刚度矩阵
% 1号区域与4号区域中a=10，2号与3号a=1
c3=((g1x.*g2x+g1y.*g2y)).*ar; % AK(1,2)=AK(2,1)=c3
c1=((g2x.*g3x+g2y.*g3y)).*ar; % AK(2,3)=AK(3,2)=c1
c2=((g3x.*g1x+g3y.*g1y)).*ar; % AK(1,3)=AK(3,1)=c2
% AK(1,1)=-AK(1,2)-AK(1,3)=-c2-c3
% AK(2,2)=-AK(2,1)-AK(2,3)=-c3-c1
% AK(3,3)=-AK(3,1)-AK(3,2)=-c1-c2
tmp=[find(t(4,:)==1), find(t(4,:)==4)];
c1(tmp)=c1(tmp)*10; c2(tmp)=c1(tmp)*10; c3(tmp)=c1(tmp)*10; 

% 装配总刚度矩阵
A=sparse(it1,it2,c3,np,np);
A=A+sparse(it2,it3,c1,np,np);
A=A+sparse(it3,it1,c2,np,np);
A=A+A.';
A=A+sparse(it1,it1,-c2-c3,np,np);
A=A+sparse(it2,it2,-c3-c1,np,np);
A=A+sparse(it3,it3,-c1-c2,np,np);

% 装配右端项
f=ar/3;
F=sparse(it1,1,f,np,1);
F=F+sparse(it2,1,f,np,1);
F=F+sparse(it3,1,f,np,1);

% 消去Dirichlet零边界
ie=zeros(1,np);
ie(e(1:2,:))=1;
ie=find(ie);
B=speye(np);
B(:,ie)=[];
A=B'*A*B;
F=B'*F;
