function [betaHeu, lambdaHeu, accHeu, betaPath, lambdaPath, accPath, tauPath, betaMeanPath, tauHeu] = sldaValHeuristic(XTrain, YTrain, XVal, YVal)
%SLDAVALHEURISTIC - heuristic validation scheme for sparse linear discriminant analysis
%sldaValHeuristic(XTrain, YTrain, XVal, YVal) performs a heuristic path validation scheme 
%for sparse linear discriminant analysis (slda) using the l1-Houdini algorithm
%
% Syntax: [betaVal, lambdaVal, accVal, plotDataVal] = sldaValFunction(XTrain, YTrain, XVal, YVal)
%
% Input:
%    XTrain - training samples from class 1 (p x nXTrain matrix)
%    YTrain - training samples from class 2 (p x nYTrain matrix)
%    XVal - validation samples from class 1 (p x nXVal matrix)
%    YVal - validation samples from class 2 (p x nYVal matrix)
%    
% Output:
%    betaHeu - linear predictor (p vector)
%    lambdaHeu - tuning parameter associated wit betaHeu (scalar)
%    accHeu - predicted accuracy (scalar)
%    betaPath - solution path of l1-HOUDINI (p x K matrix)
%    lambdaPath - associated values of the homotopy parameter (K vector)
%    accPath - associated accuracy values w.r.t. validation data (K vector)
%    tauPath - associated optimal thresholds (K vector)
%    betaMeanPath - associated values of beta' * mu (K vector)
%    tauHeu - best accuracy on the path (scalar)
%
% Example: see exampleFSC.m
%
% The vector betaHeu can be used to classify samples z with unknown class
% via a variant of Fisher's linear discriminant rule. Writing
% mu := (mean(XTrain, 2) + mean(YTrain, 2))/2, the standard way to classify
% a new sample would be to assign it to class 1 if and only if
% z' * betaHeu >= mu' * betaHeu. This function selects a threshold tau
% such that the classification criterion z' * betaHeu >= tau is optimal
% w.r.t. the validation data.
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,
% dualActiveSet.m, sldaThreshold.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET, SLDATHRESHOLD

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 04-March-2019

% use l1-HOUDINI to compute solution path on training data
XBar = mean(XTrain, 2);
YBar = mean(YTrain, 2);
SigmaHat = ((XTrain - XBar) * (XTrain - XBar)' + ...u
    (YTrain - YBar) * (YTrain - YBar)') / (size(XTrain, 2) + size(YTrain, 2));
[betaPath, lambdaPath] = l1HoudiniPath(SigmaHat, XBar - YBar, 0);
K = length(lambdaPath);

% compute accuracies along solution path with respective best thresholds
tauPath = zeros(K, 1);
accPath = zeros(K, 1);
for k = 1:K
    [~, ~, tauPath(k), accPath(k)] = sldaThreshold(betaPath(:, k), XVal, YVal);
end

[accHeu, indexHeu] = max(accPath);
lambdaHeu = lambdaPath(indexHeu);
betaHeu = betaPath(:, indexHeu);
tauHeu = tauPath(indexHeu);
betaMeanPath = betaPath' * (XBar + YBar) / 2;