function [x,y]=square0101g(bs,s)

nbs=4;

if nargin==0
  x=nbs; % number of boundary segments
  return
end

d=[
  0 0 0 0 % start parameter value
  1 1 1 1 % end parameter value
  0 0 0 0 % left hand region
  1 1 1 1 % right hand region
];

bs1=bs(:)';

if find(bs1<1 | bs1>nbs)
  error(message('pde:squareg:InvalidBs'))
end

if nargin==1
  x=d(:,bs1);
  return
end

x=zeros(size(s));
y=zeros(size(s));
[m,n]=size(bs);
if m==1 && n==1
  bs=bs*ones(size(s)); % expand bs
elseif m~=size(s,1) || n~=size(s,2)
  error(message('pde:squareg:SizeBs'));
end

if ~isempty(s)

% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=interp1([d(1,1),d(2,1)],[0 1],s(ii));
y(ii)=interp1([d(1,1),d(2,1)],[1 1],s(ii));
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=interp1([d(1,2),d(2,2)],[1 1],s(ii));
y(ii)=interp1([d(1,2),d(2,2)],[1 0],s(ii));
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=interp1([d(1,3),d(2,3)],[1 0],s(ii));
y(ii)=interp1([d(1,3),d(2,3)],[0 0],s(ii));
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=interp1([d(1,4),d(2,4)],[0 0],s(ii));
y(ii)=interp1([d(1,4),d(2,4)],[0 1],s(ii));
end

end

