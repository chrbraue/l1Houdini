function [betaCV, lambdaCV, accCV, plotDataCV] = sldaCVFunction(X, Y, N)
%SLDACVFUNCTION - cross validation scheme for sparse linear discriminant analysis
%sldaCVFunction(X, Y, N) performs an N-fold cross-validation scheme for
%sparse linear discriminant analysis (slda) using the l1-Houdini algorithm
%
% Syntax: [betaCV, lambdaCV, accCV, plotDataCV] = sldaCVFunction(X, Y, N)
%
% Input:
%    X - samples from class 1 (p x nX matrix)
%    Y - samples from class 2 (p x nY matrix)
%    N - number of folds (integer)
%    
% Output:
%    betaCV - linear predictor (p vector)
%    lambdaCV - tuning parameter associated wit betaCV (scalar)
%    accCV - predicted accuracy (scalar)
%    plotDataCV - cell containing plot data (estimated out-of-sample accuracy)
%
% Example: see exampleFSC.m
%
% The vector betaCV can be used to classify samples z with unknown class
% via a variant of Fisher's linear discriminant rule. Writing
% mu := (mean(X, 2) + mean(Y, 2))/2, the new sample is assigned to class 1
% if and only if (z - mu)' * betaCV >= 0.
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,  dualActiveSet.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 04-March-2019
    
% get numbers of positive / negative training examples
nX = size(X, 2);
nY = size(Y, 2);    
% error if #(training examples per class) < #(folds)
assert(nX >= N, 'number of training examples in class 1 is smaller than number of folds');
assert(nY >= N, 'number of training examples in class 2 is smaller than number of folds');
% calculate fold sizes (only last folds may be smaller)
foldSizeX = floor(nX / N);
foldSizeY = floor(nY / N);    
% generate random folds
permX = randperm(nX);
permY = randperm(nY);    
% allocate memory
muHat = cell(1, N);
betaPathsHoudini = cell(1, N);
lambdaPathsHoudini = cell(1, N);
lambdaPathsCV = cell(1, N);
lambdaPathsCVUnion = [];
accuracyPathsCV = cell(1, N);
lambdaMin = 0;
lambdaMax = 0;

% calculate each folds solution and parameter paths
for l = 1:N
    % determine training data for current fold
    if l < N
        trainFoldX = permX([1:(l - 1) * foldSizeX, l * foldSizeX + 1:end]);
        trainFoldY = permY([1:(l - 1) * foldSizeY, l * foldSizeY + 1:end]);
    else
        trainFoldX = permX(1:(l - 1) * foldSizeX);
        trainFoldY = permY(1:(l - 1) * foldSizeY);
    end
    % calculate empirical means and covariance matrix
    XBar = mean(X(:, trainFoldX), 2);
    YBar = mean(Y(:, trainFoldY), 2);
    muHat{l} = (XBar + YBar) / 2;
    SigmaHat = ((X(:, trainFoldX) - XBar) * (X(:, trainFoldX) - XBar)'  ...
        + (Y(:, trainFoldY) - YBar) * (Y(:, trainFoldY) - YBar)') ...
        / (length(trainFoldX) + length(trainFoldY));
    % calculate l-th solution path
    [betaPathsHoudini{l}, lambdaPathsHoudini{l}] = ...
        l1HoudiniPath(SigmaHat, XBar - YBar, lambdaMin);
    lambdaMin = max(lambdaMin, lambdaPathsHoudini{l}(end));
    lambdaMax = max(lambdaMax, lambdaPathsHoudini{l}(1));
end
    
% calculate each folds parameter and accuracy paths
for l = 1:N
    % determine testing data for current fold
    if l < N
        testFoldX = permX((l - 1) * foldSizeX + 1:l * foldSizeX);
        testFoldY = permY((l - 1) * foldSizeY + 1:l * foldSizeY);
    else
        testFoldX = permX((l - 1) * foldSizeX + 1:end);
        testFoldY = permY((l - 1) * foldSizeY + 1:end);
    end
    % initialize paths (if beta = 0, exactly all samples from class 1 are
    % classified correctly)
    lambdaPathCV = lambdaPathsHoudini{l}(1);
    accuracyPathCV = length(testFoldX);
    % compute jumps of CV_l along the path up to lambdaMin
    for k = 2:length(lambdaPathsHoudini{l})
        % calculate direction
        d = betaPathsHoudini{l}(:, k) - betaPathsHoudini{l}(:, k - 1);            
        % calculate jump discontinuities w.r.t. X
        g = (X(:, testFoldX) - muHat{l})' * d;
        s = -(X(:, testFoldX) - muHat{l})' * ...
            betaPathsHoudini{l}(:, k - 1) ./ g;
        cX = (g > 0 & s > 0 & s < 1) - (g < 0 & s >= 0 & s < 1);
        lambdaX = lambdaPathsHoudini{l}(k - 1) + s(cX ~= 0) * ...
            (lambdaPathsHoudini{l}(k) - lambdaPathsHoudini{l}(k - 1));
        cX = cX(cX ~= 0);
        % calculate jump discontinuities w.r.t. Y
        h = (Y(:, testFoldY) - muHat{l})' * d;
        t = -(Y(:, testFoldY) - muHat{l})' * ...
            betaPathsHoudini{l}(:, k - 1) ./ h;
        cY = (h < 0 & t >= 0 & t < 1) - (h > 0 & t > 0 & t  < 1);
        lambdaY = lambdaPathsHoudini{l}(k - 1) + t(cY ~= 0) * ...
            (lambdaPathsHoudini{l}(k) - lambdaPathsHoudini{l}(k - 1));
        cY = cY(cY ~= 0);
        % compute union of jumps w.r.t. X and Y
        [lambdaTmp, indexTmp] = sort([lambdaX; lambdaY]);
        lambdaPathCV = [lambdaTmp; lambdaPathCV]; %#ok<*AGROW>
        cTmp = [cX; cY];
        accuracyPathCV = [accuracyPathCV(1) + cumsum(cTmp(indexTmp), 'reverse')...
            ; accuracyPathCV];
    end    
    % remove redundant entries and add one for lambdaMax + eps (so that
    % accuracy for lambda > lambdaMax is represented)
    [lambdaPathsCV{l}, indexTmp] = unique(lambdaPathCV);
    lambdaPathsCV{l} = [lambdaPathsCV{l}; lambdaMax + (lambdaMax - lambdaMin) / 10];
    accuracyPathsCV{l} = [accuracyPathCV(indexTmp); length(testFoldX)] / ...
        (length(testFoldX) + length(testFoldY));
    % add entry for lambdaMin (which will be the smallest point on the
    % final path) and remove entries below lambdaMin
    [lambdaPathsCVTmp, ~, indexTmp] = unique([lambdaMin; lambdaPathsCV{l}]);
    if indexTmp(1) > 1
        accuracyPathsCV{l} = [accuracyPathsCV{l}(1:indexTmp(1) - 1); ...
            interp1(lambdaPathsCV{l}, accuracyPathsCV{l}, lambdaMin, 'next'); ...
            accuracyPathsCV{l}(indexTmp(1):end)];
        lambdaPathsCV{l} = lambdaPathsCVTmp;
        indexTmp2 = lambdaPathsCV{l} >= lambdaMin;
        lambdaPathsCV{l} = lambdaPathsCV{l}(indexTmp2);
        accuracyPathsCV{l} = accuracyPathsCV{l}(indexTmp2);
    else
        lambdaPathsCV{l} = [lambdaMin; lambdaPathsCV{l}];
        accuracyPathsCV{l} = [accuracyPathsCV{l}(1); accuracyPathsCV{l}];
    end
    % add l-th parameter path to final parameter grid
    lambdaPathsCVUnion = unique([lambdaPathsCVUnion; lambdaPathsCV{l}]);    
end
    
% calculate mean of all accuracy curves
accuracyPathsInterp = zeros(length(lambdaPathsCVUnion), N);
for l = 1:N
    accuracyPathsInterp(:, l) = interp1(lambdaPathsCV{l}, accuracyPathsCV{l}, ...
        lambdaPathsCVUnion, 'next');
end
accuracyPathsMean = mean(accuracyPathsInterp, 2);
    
% calculate optimal lambda
[accCV, indexAccuracyCV] = max(accuracyPathsMean);
lambdaCV = lambdaPathsCVUnion(indexAccuracyCV);
if lambdaCV > lambdaMin
   lambdaCV = (lambdaCV + lambdaPathsCVUnion(indexAccuracyCV - 1)) / 2; 
end
    
% calculate final model with best parameter and full data
XBar = mean(X, 2);
YBar = mean(Y, 2);
SigmaHat = ((X - XBar) * (X - XBar)' + (Y - YBar) * (Y - YBar)') / (nX + nY);
betaTmp = l1HoudiniPath(SigmaHat, XBar - YBar, lambdaCV);
betaCV = betaTmp(:, end);
    
% prepare struct with plot data
plotDataCV = {};
plotDataCV{1} = [lambdaMin; lambdaPathsCVUnion];
plotDataCV{2} = [accuracyPathsMean; size(X, 2) / (size(X, 2) + size(Y, 2))];
plotDataCV{3} = lambdaMin;
plotDataCV{4} = lambdaMax;
plotDataCV{5} = accuracyPathsInterp;