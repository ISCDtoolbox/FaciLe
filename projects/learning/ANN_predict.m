function [y] = ANN_predict(aNN,Xtest)

if aNN.flag_standardize==1
    Xtest = bsxfun(@minus,Xtest, aNN.mu);
    Xtest = bsxfun(@rdivide,Xtest,aNN.we);
end
hatytest  = aNN.net(Xtest');
y = hatytest';
%hatytest = vec2ind(hatytest)'-1;
end