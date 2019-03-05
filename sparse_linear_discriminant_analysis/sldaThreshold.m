function [tau, acc, tauBest, accBest] = sldaThreshold(beta, XVal, YVal)
%SLDATHRESHOLD - auxiliary function for optimal threshold selection
%sldaThreshold(beta, XVal, YVal) performs an efficient enumeration scheme
%to select a threshold tau optimally w.r.t to validation data when the
%classification criterion z' * beta >= tau is to be applied.
%
% Syntax: [tau, acc, tauBest, accBest] = sldaThreshold(beta, XVal, YVal)
%
% Input:
%    beta - linear predictor (p vector)
%    XVal - validation samples from class 1 (p x nXVal matrix)
%    YVal - validation samples from class 2 (p x nYVal matrix)
%    
% Output:
%    tau - vector of tested tau values
%    acc - associated accuracies on the validation data
%    tauBest - best tau
%    accBest - best accuracy
%
% Example: see exampleFSC.m
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: SLDAVALHEURISTIC

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 04-March-2019

betaX = sort(XVal' * beta); nX = length(betaX);
betaY = sort(YVal' * beta); nY = length(betaY);

tau = betaX(1);
[~, indexYStart] = min(betaY < tau);
indexYStart = indexYStart - 1;
acc = (nX + indexYStart);

indexX = 1;
for indexY = indexYStart + 1:nY
    tau = [tau; betaY(indexY) + 1e-10]; %#ok<*AGROW>
    accTmp = acc(end) + 1;
    while indexX <= nX && betaX(indexX) < tau(end)
        indexX = indexX + 1;
        accTmp = accTmp - 1;
    end
    acc = [acc; accTmp];
end

acc = acc / (nX + nY);

[accBest, indexBest] = max(acc);
tauBest = tau(indexBest);