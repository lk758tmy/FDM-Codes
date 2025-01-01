function I=mg_prolangation(p1,p2,t1)
% 构造延拓矩阵
% 由于refinemesh是将新的节点追加到原有的后面，构造延拓矩阵是容易的

tSize=size(t1,2);
height=size(p1,2);
width=size(p2,2);

% 本机为4核，手动实现并行化
I4={};
p={t1(1:3,1:(tSize/4)) t1(1:3,(tSize/4+1):(2*tSize/4)) ...
   t1(1:3,(2*tSize/4+1):(3*tSize/4)) t1(1:3,(3*tSize/4+1):tSize)};
parfor j=1:4
    I4{j}=mg_prolangation_sub(height,width,p{j});
end
I=max(max(I4{1},I4{2}),max(I4{3},I4{4}));
I=[speye(width);I];