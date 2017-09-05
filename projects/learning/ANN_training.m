function [aNN]= ANN_training(Xtrain,ytrain,flag_standardize)
%builder for a neural network for regression
%Xtrain= ntrain x P matrix
%ytrain= ntrain x 1 -> it must contain only 0s and 1s


%%
%step 0: initialization
if nargin==2
    flag_standardize=0;
end
    
%flag_standardize=0;

trainFcn = 'trainlm';  %algorithm used for training, see fitnet for more options
%trainFcn ='trainbr';

hiddenLayerSize = 10; %number of neurons in the first hidden layer (if one wants more layers then use [10,5,2,...])

trainRatio=80/100;
valRatio = 20/100;
testRatio = 0/100;

if flag_standardize==1
    [Ztrain, aNN.mu, aNN.we] = zscore(Xtrain);
    aNN.we(aNN.we==0) = 1;
else
    Ztrain=Xtrain;
end
Ztrain=Ztrain';
%ytrain=dummyvar(ytrain+1)'; %this is for classification

%%
% Create a Pattern Recognition Network
net = fitnet(hiddenLayerSize,trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio =trainRatio; 
net.divideParam.valRatio =  valRatio;
net.divideParam.testRatio = testRatio;

%%
% Train the Network
net = train(net,Ztrain,ytrain');
aNN.net=net;
aNN.flag_standardize=flag_standardize;
end