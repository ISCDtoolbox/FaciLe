function [Xlearn,ylearn,Xtest,ytest]= random_partition(X,y,ntrain_rel)
%this function partitions the dataset randomly
%ntrain_rel = percentage of data in the learning set

%%
N=size(X,1);
ntrain=round(N*ntrain_rel);
I=randperm(N);

%%
Xlearn=X(I(1:ntrain),:);
ylearn=y(I(1:ntrain),:);
Xtest=X(I(ntrain+1:end),:);
ytest=y(I(ntrain+1:end),:);


end