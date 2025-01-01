%%%
% Problem:
%     -Gradient(a*Gradient(u))=1 in Omega=(-1,1)x(-1,1)
%     a=1 x*y>0 / a=10 x*y<=0
%     u(x,y)=0 on the border of Omega
%%%

m=10; %磨光次数
% l=20; %完全多重网格每层迭代次数
% tol=1e-6; %误差限

%生成区域和初始网格
gd=[[3	3	3	3]
    [4	4	4	4]
    [0	-1	-1	0]
    [1	0	0	1]
    [1	0	0	1]
    [0	-1	-1	0]
    [1	1	0	0]
    [1	1	0	0]
    [0	0	-1	-1]
    [0	0	-1	-1]];
% 用四个正方形拼起来，保证a(x,y)的不连续界面是单元的边
% 但在initmesh时会生成"多余"的区域边界(即两个正方形的公共边)
g=decsg(gd);
[p,e,t]=initmesh(g,"Hmax",1); %很粗的网格

% 删去"多余"的区域边界
tmp=[];
for j=1:size(e,2)
    if(((abs(abs(p(1,e(1,j)))-1)>1e-8)&&(abs(abs(p(2,e(1,j)))-1)>1e-8))|| ...
       ((abs(abs(p(1,e(2,j)))-1)>1e-8)&&(abs(abs(p(2,e(2,j)))-1)>1e-8)))
        tmp=[tmp,j];
    end
end
e(:,tmp)=[];

% 装配各级刚度矩阵和右端向量
A={}; F={}; B={};
p={p}; e={e}; t={t};
for i=1:8
    if i>1
        [p{i},e{i},t{i}]=refinemesh(g,p{i-1},e{i-1},t{i-1});
    end
    % pdemesh(p,e,t)
    [A{i},F{i},B{i}]=pdeasmpoi(p{i},e{i},t{i});
end

% 构造延拓矩阵
I={};
for i=2:8
    fprintf("构造第%d层延拓矩阵",i-1)
    tic;
    I{i-1}=B{i}'*mg_prolangation(p{i},p{i-1},t{i})*B{i-1};
    toc
end

fprintf("\n")

% 求解pde
tol=[0 0 0];
for i=2:8
    fprintf("第%d加密，内部结点个数：%d\n",i-1,size(A{i},1));

    fprintf("--直接求解--\n");
    tic;
    % uDirect=full(B{i}*(A{i}\F{i}));
    % pdesurf(p{i},t{i},uDirect)
    uDirect=A{i}\F{i};
    toc

    fprintf("--完全网格法--\n");
    for l=1:6
        fprintf("l=%d\n",l*5);
        tic;
        uFull=fmg(A,F,B,I,l*5,m,i);
        toc
        fprintf("与直接法的差距(L∞)：%e\n",0+max(abs(uDirect-uFull)))
        % 此处的+0是为了避免"没有为稀疏输入定义函数"报错
        if(mod(l,2)==0)
            tol(l/2)=(max(abs(F{i}-A{i}*uFull))+0)/(max(abs(F{i}))+0);
            %为了可比性，V循环的残差限使用完全网格法最终的残差
        end
    end

    fprintf("--V循环算法--\n");
    for l=1:3
        fprintf("tol=%.1e (对应完全网格法l=%d)\n",tol(l),l*10);
        tic;
        [uV,stepsV]=mg_vcycle(A{i},F{i},I,sparse(size(A{i},1),1),m,i,tol(l));
        toc
        fprintf("迭代次数：%d\n",stepsV)
        fprintf("与直接法的差距(L∞)：%e\n",0+max(abs(uDirect-uV)))
    end

    fprintf("\n")
end
