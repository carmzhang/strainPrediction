function [out,f] = svmdecision(Xnew,svm_struct)
%SVMDECISION Evaluates the SVM decision function

%   Copyright 2004-2012 The MathWorks, Inc.

for c = 1:size(Xnew, 2)
    Xnew(:,c) = svm_struct.ScaleData.scaleFactor(c) * ...
    (Xnew(:,c) +  svm_struct.ScaleData.shift(c));
end
sv = svm_struct.SupportVectors;
alphaHat = svm_struct.Alpha;
bias = svm_struct.Bias;
kfun = svm_struct.KernelFunction;
kfunargs = svm_struct.KernelFunctionArgs;

f = (feval(kfun,sv,Xnew,kfunargs{:})'*alphaHat(:)) + bias;
out = sign(f);
% points on the boundary are assigned to class 1
out(out==0) = 1;
