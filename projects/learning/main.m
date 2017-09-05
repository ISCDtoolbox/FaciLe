%16/05/2017
%main file for multiple outputs regression
%below I use neural networks

%% I load the data 
%X(i,:)= [x(i),y(i),\mu(i)]
%Y(i,:) = [real(u(x,y,\mu)), imag(u(x,y,\mu))]
% load('data_helmhotz.mat')
% nlearn = 0.1; %percentage of points used for learning (train+validation)
% flag_standardize=1; %it is common practice to standardize input data (it never hurts, 
%                     %and sometimes it helps)
%                     %it might be worth to standardize the output (since we
%                     %have more than one output), but for now I don't do it
% d = size(data.Y,2);

% Load Data Masseter
dataSet = cell(1, 2);
dataSet{2} = alphaNoRot;
dataSet{1} = dlmread('Mesures.csv', ';');
dataSet{1} = dataSet{1}.';
dataSet{2}(:,17)= []; % Pas de mesures associé à ces coefficient de masséter car perte du CTscan
dataSet{2}(:,18)= [];
for i = 1:64
    dataSet{2}(1,:) = [];
end
% 
% nlearn = 0.9;
% flag_standardize=1;
% d = size(dataSet{1},2);


%% Ma version simplifiée
net = fitnet(50);
%view(net)
[net, tr] = train(net, dataSet{1}, dataSet{2});
nntraintool
a = net(dataSet{1}(:,1));
%
coef = zeros(5,68);
for i = 1:68
    coef(:,i) = net(dataSet{1}(:,i));
end

%% Reconstruction masséter
n = 5;
numFiles = 10;

RRecoNoRot = cell(n,numFiles);
RecoNoRot = cell(1, numFiles);
rmeshNoRot = cell(1, numFiles);
rsolNoRot = cell(1, numFiles);

for i = 1:numFiles
    RecoNoRot{i} =  zeros(numLines,3);
    rmeshNoRot{i} = zeros(numLines,3);
    rsolNoRot{i} = zeros(numLines,3);
end

for i = 1:n %Numero des vecteur pris en compte
    for j = 1:numFiles % Pour tous les masseters + le 52
        RRecoNoRot{i,j} = coef(i,j)*BBNoRot{i+64};
        RecoNoRot{j} = RecoNoRot{j} + RRecoNoRot{i,j};
        rmeshNoRot{j} = template + RecoNoRot{j};
        rsolNoRot{j} = LNoRot{j} - RecoNoRot{j};
    end
end

for fileNum = 1:numFiles
    fileName = sprintf('Reco%dNoRot.%d.mesh',n,fileNum);
    okkmeshNoRot = writemesh(fileName, rmeshNoRot{fileNum}.', tri,edg,crn);
    fileName = sprintf('Reco%dNoRot.%d.sol',n,fileNum);
    okksolNoRot = writesol(fileName, rsolNoRot{fileNum}.');
end

%% I perform the training
%the function ANN_training depends on the number of the number of layers
%used. To be fair, I have no experience with how to set this number. The default is 10
%but the general practice is to increase the number of layers and neurons
%per layer as the number of training points increases.

% [Xlearn,ylearn,Xtest,ytest]= random_partition(data.X,data.Y,nlearn);
% ntest = size(Xtest,1);
% fprintf('   --> neural nets ...')
% [aNN]= ANN_training(Xlearn,ylearn,flag_standardize);

%Moi
[Xlearn,ylearn,Xtest,ytest]= random_partition(dataSet{1},dataSet{2},nlearn);
ntest = size(Xtest,1);
fprintf('   --> neural nets ...')
[aNN]= ANN_training(Xlearn,ylearn,flag_standardize);

%% I assess performance by computing the out-of-sample R^2 (coefficient of determination)
%(this is a reasonable indicator, but it is definitely not the only one)
barytrain = mean(ylearn,1);
haty = ANN_predict(aNN,Xtest);
R2 = 1 - mean( sum((haty - ytest).^2,2) )/mean( sum( bsxfun(@minus, ytest, barytrain).^2,2) )

fprintf(' Done.\n\n');


%% below, I provide results for another method (CART, or regression trees)
%here I train separately the two components. Performance seems comparable
%for this dataset.
fprintf('   --> CART ...')
%CARTs (at least the Matlab implementation) do not work with multiple
%outputs
for n=1:d
    rT{n} = fitrtree(Xlearn,ylearn(:,n));
end

haty= zeros( ntest,d);
for n=1:d
    haty(:,n) = predict(rT{n},Xtest);
end

R2 = 1 - mean( sum((haty - ytest).^2,2) )/mean( sum( bsxfun(@minus, ytest, barytrain).^2,2) )

fprintf(' Done.\n\n');
% 
% 
% 



