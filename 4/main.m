%%%
% Problem:
%     Laplasian(u)=-1 in (0,1)x(0,1)
%     u(x,0)=u(x,1)=0
%     u(0,y)=u(1,y)
%%%

g=@square0101g;
[p,e,t]=poimesh(g);
[p,e,t]=refinemesh(g,p,e,t);
[p,e,t]=refinemesh(g,p,e,t); % start with a=0.25 (32 triangles)
% pdemesh(p,e,t);

uexact=(0.5*p(2,:).*(1-p(2,:)))';
u=pdeasmpoi_period(p,e,t);
% pdesurf(p,t,u-uexact);

size(p(1,:))
max(abs(u-uexact))
sqrt(sum(0.25*(u-uexact).^2))
for i=1:5
    [p,e,t]=refinemesh(g,p,e,t);
    uexact=(0.5*p(2,:).*(1-p(2,:)))';
    u=pdeasmpoi_period(p,e,t);
    size(p(1,:))
    max(abs(u-uexact))
    sqrt(sum(2^(-i-2)*(u-uexact).^2))
end