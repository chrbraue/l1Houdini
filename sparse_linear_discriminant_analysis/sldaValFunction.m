function [betaVal, lambdaVal, accVal, plotDataVal] = sldaValFunction(XTrain, YTrain, XVal, YVal)
%SLDAVALFUNCTION - validation scheme for sparse linear discriminant analysis
%sldaValFunction(XTrain, YTrain, XVal, YVal) performs a path validation scheme for 
%sparse linear discriminant analysis (slda) using the l1-Houdini algorithm
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
%    betaVal - linear predictor (p vector)
%    lambdaVal - tuning parameter associated wit betaVal (scalar)
%    accVal - predicted accuracy (scalar)
%    plotDataVal - cell containing plot data (estimated out-of-sample accuracy)
%
% Example: see exampleFSC.m
%
% The vector betaVal can be used to classify samples z with unknown class
% via a variant of Fisher's linear discriminant rule. Writing
% mu := (mean(XTrain, 2) + mean(YTrain, 2))/2, the new sample is assigned
% to class 1 if and only if (z - mu)' * betaVal >= 0.
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,  dualActiveSet.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 04-March-2019

% use l1-HOUDINI to compute solution path on training data
XBar = mean(XTrain, 2);
YBar = mean(YTrain, 2);
SigmaHat = ((XTrain - XBar) * (XTrain - XBar)' + ...
    (YTrain - YBar) * (YTrain - YBar)') / (size(XTrain, 2) + size(YTrain, 2));
[betaPathHoudini, lambdaPathHoudini] = l1HoudiniPath(SigmaHat, XBar - YBar, 0);

% get numbers of positive / negative training examples in validation data
nXVal = size(XVal, 2);
nYVal = size(YVal, 2);    
lambdaMin = lambdaPathHoudini(end);
lambdaMax = lambdaPathHoudini(1);

% initialize paramter path and accuracy path
lambdaPathVal = lambdaPathHoudini(1);
accuracyPathVal = nXVal;

% compute jumps of accuracy path
muHat = (XBar + YBar) / 2;
for k = 2:length(lambdaPathHoudini)
    % calculate direction
    d = betaPathHoudini(:, k) - betaPathHoudini(:, k - 1);            
    % calculate jump discontinuities w.r.t. X
    g = (XVal - muHat)' * d;
    s = -(XVal - muHat)' * betaPathHoudini(:, k - 1) ./ g;
    cX = (g > 0 & s > 0 & s < 1) - (g < 0 & s >= 0 & s < 1);
    lambdaX = lambdaPathHoudini(k - 1) + s(cX ~= 0) * ...
        (lambdaPathHoudini(k) - lambdaPathHoudini(k - 1));
    cX = cX(cX ~= 0);
    % calculate jump discontinuities w.r.t. Y
    h = (YVal - muHat)' * d;
    t = -(YVal - muHat)' * betaPathHoudini(:, k - 1) ./ h;
    cY = (h < 0 & t >= 0 & t < 1) - (h > 0 & t > 0 & t  < 1);
    lambdaY = lambdaPathHoudini(k - 1) + t(cY ~= 0) * ...
        (lambdaPathHoudini(k) - lambdaPathHoudini(k - 1));
    cY = cY(cY ~= 0);
    % compute union of jumps w.r.t. X and Y
    [lambdaTmp, indexTmp] = sort([lambdaX; lambdaY]);
    lambdaPathVal = [lambdaTmp; lambdaPathVal]; %#ok<*AGROW>
    cTmp = [cX; cY];
    accuracyPathVal = [accuracyPathVal(1) + cumsum(cTmp(indexTmp), 'reverse');
        accuracyPathVal];
end

% remove redundant entries and add one for values greater than lambdaMax
[lambdaPathVal, indexTmp] = unique(lambdaPathVal);
lambdaPathVal = [lambdaPathVal; lambdaMax + (lambdaMax - lambdaMin) / 10];
accuracyPathVal = [accuracyPathVal(indexTmp); nXVal] / (nXVal + nYVal);

% add entry for lambdaMin
if lambdaPathVal(1) > lambdaMin
    lambdaPathVal = [lambdaMin; lambdaPathVal];
    accuracyPathVal = [accuracyPathVal(1); accuracyPathVal];
end
    
% calculate optimal lambda and best predictor
[accVal, indexAccuracyVal] = max(accuracyPathVal);
lambdaVal = lambdaPathVal(indexAccuracyVal);
if lambdaVal > lambdaMin
    lambdaVal = (lambdaVal + lambdaPathVal(indexAccuracyVal - 1)) / 2;
end
betaVal = l1Houdini(SigmaHat, XBar - YBar, lambdaVal);
  
% prepare structs with plot data
plotDataVal = {};
plotDataVal{1} = [lambdaMin; lambdaPathVal];
plotDataVal{2} = [accuracyPathVal; nXVal / (nXVal + nYVal)];
plotDataVal{3} = lambdaMin;
plotDataVal{4} = lambdaMax;