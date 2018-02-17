%SLDACVEXAMPLE - random example for sldaCVFunction
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,
%       dualActiveSet.m, sldaCVFunction.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET,  SLDACVFUNCTION

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 17-February-2018

% set dimension of feature space
p = 100;

% set number of samples from classes 1 and 2
nX = 75;
nY = 75;

% generate random positive definite covariance matrix
tmp = randn(p); Sigma = eye(p) + tmp' * tmp;

% set means of distributions 1 and 2
muX = ones(p, 1);
muY = -ones(p, 1);

% draw random samples from distributions 1 and 2
X = mvnrnd(muX, Sigma, nX)';
Y = mvnrnd(muY, Sigma, nY)';

% run cross validation scheme
[betaCV, lambdaCV, rateCV, figureHandle] = sldaCVFunction(X, Y, 10, true);

% comute success rate training data
XBar = mean(X, 2);
YBar = mean(Y, 2);
muHat = (XBar + YBar) / 2;
rateTrain = (sum((X - muHat)' * betaCV >= 0) + ...
    sum((Y - muHat)' * betaCV < 0)) / (nX + nY);

% compute success rate on independent testing data
XTest = mvnrnd(muX, Sigma, nX)';
YTest = mvnrnd(muY, Sigma, nY)';
rateTest = (sum((XTest - muHat)' * betaCV >= 0) + ...
    sum((YTest - muHat)' * betaCV < 0)) / (nX + nY);

% display results
fprintf('\nselected parameter:\t\t lambdaCV = %f\n', lambdaCV);
fprintf('predicted success rate:\t\t rateCV = %f\n', rateCV);
fprintf('success rate on training set:\t rateTrain = %f\n', rateTrain);
fprintf('success rate on testing set:\t rateTest = %f\n', rateTest);
set(figureHandle, 'Visible', 'on');