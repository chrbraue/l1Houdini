function [betaCV, lambdaCV, rateCV, figureHandle] = sldaCV(X, Y, N, plot)
%SLDACV - cross validation scheme for sparse linear discriminant analysis
%sldaCV(X, Y, N, plot) performs an N-fold cross-validation scheme for
%sparse linear discriminant analysis (slda) using the l1-Houdini homotopy
%algorithm.
%
% Syntax: [betaCV, lambdaCV, rateCV, figureHandle] = sldaCV(X, Y, N, plot)
%
% Input:
%    X - samples from class 1 (p x nX matrix)
%    Y - samples from class 2 (p x nY matrix)
%    N - number of folds (integer)
%    plot - switches plot on/off (boolean)
%    
% Output:
%    betaCV - linear predictor (p vector)
%    lambdaCV - tuning parameter associated wit betaCV
%    rateCV - predicted success rate
%    figureHandle - figure handle (if plot = true) or empty (else)
%
% Example:
% p = 100; nX = 75; nY = 75; (data dimension, numbers of samples)
% tmp = randn(p); Sigma = eye(p) + tmp' * tmp; (pos. def. covariance mat.)
% mu1 = ones(p, 1); mu2 = -ones(p, 1); (means of distributions)
% X = mvnrnd(mu1, Sigma, nX)'; Y = mvnrnd(mu2, Sigma, nY)'; (data)
% [betaCV, lambdaCV, rateCV, figureHandle] = sldaCV(X, Y, 10, true);
% set(figureHandle, 'Visible', 'on');
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
% February 2018; Last revision: 15-February-2018
    
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
ratePathsCV = cell(1, N);
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
        l1HoudiniPath(SigmaHat, XBar - YBar, 0);
    lambdaMin = max(lambdaMin, lambdaPathsHoudini{l}(end));
    lambdaMax = max(lambdaMax, lambdaPathsHoudini{l}(1));
end
    
% calculate each folds parameter and success rate paths
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
    ratePathCV = length(testFoldX);
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
        ratePathCV = [ratePathCV(1) + cumsum(cTmp(indexTmp), 'reverse')...
            ; ratePathCV];
    end    
    % remove redundant entries and add one for lambdaMax + eps (so that
    % success rate for lambda > lambdaMax is represented)
    [lambdaPathsCV{l}, indexTmp] = unique(lambdaPathCV);
    lambdaPathsCV{l} = [lambdaPathsCV{l}; lambdaMax + 1];
    ratePathsCV{l} = [ratePathCV(indexTmp); length(testFoldX)] / ...
        (length(testFoldX) + length(testFoldY));
    % add entry for lambdaMin (which will be the smallest point on the
    % final path) and remove entries below lambdaMin
    [lambdaPathsCVTmp, ~, indexTmp] = unique([lambdaMin; lambdaPathsCV{l}]);
    if indexTmp(1) > 1
        ratePathsCV{l} = [ratePathsCV{l}(1:indexTmp(1) - 1); ...
            interp1(lambdaPathsCV{l}, ratePathsCV{l}, lambdaMin, 'next'); ...
            ratePathsCV{l}(indexTmp(1):end)];
        lambdaPathsCV{l} = lambdaPathsCVTmp;
        indexTmp2 = lambdaPathsCV{l} >= lambdaMin;
        lambdaPathsCV{l} = lambdaPathsCV{l}(indexTmp2);
        ratePathsCV{l} = ratePathsCV{l}(indexTmp2);
    else
        lambdaPathsCV{l} = [lambdaMin; lambdaPathsCV{l}];
        ratePathsCV{l} = [ratePathsCV{l}(1); ratePathsCV{l}];
    end
    % add l-th parameter path to final parameter grid
    lambdaPathsCVUnion = unique([lambdaPathsCVUnion; lambdaPathsCV{l}]);    
end
    
% calculate mean of all success rate curves
ratePathsInterp = zeros(length(lambdaPathsCVUnion), N);
for l = 1:N
    ratePathsInterp(:, l) = interp1(lambdaPathsCV{l}, ratePathsCV{l}, ...
        lambdaPathsCVUnion, 'next');
end
ratePathsMean = mean(ratePathsInterp, 2);
    
% calculate optimal lambda
[rateCV, indexRateCV] = max(ratePathsMean);
lambdaCV = lambdaPathsCVUnion(indexRateCV);
if lambdaCV > lambdaMin
   lambdaCV = (lambdaCV + lambdaPathsCVUnion(indexRateCV - 1)) / 2; 
end
    
% calculate final model with best parameter and full data
XBar = mean(X, 2);
YBar = mean(Y, 2);
SigmaHat = ((X - XBar) * (X - XBar)'  ...
        + (Y - YBar) * (Y - YBar)') ...
        / (nX + nY);
betaTmp = l1HoudiniPath(SigmaHat, XBar - YBar, lambdaCV);
betaCV = betaTmp(:, end);
    
% plot ratePathsMean
figureHandle = [];
if plot
    figureHandle = figure;
    set(gcf, 'Visible', 'off');
    set(gca, 'XLim', [lambdaMin, lambdaMax + 1]);
    set(gca, 'YLim', [min(min(ratePathsMean), .45), 1]);
    xlabel('$\lambda$', 'Interpreter', 'latex');
    ylabel('$\mathrm{CV}(\lambda)$', 'Interpreter', 'latex');
    hold on;
    stairs([lambdaMin; lambdaPathsCVUnion], [ratePathsMean; 0.5], 'k');
end
