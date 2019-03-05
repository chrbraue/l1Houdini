%EXAMPLEMNISTSPARSITY - reproduces Fig. 6 in (*)
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,
%       dualActiveSet.m, sldaValFunction.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET,  SLDAVALFUNCTION
%
% (*) C. Brauer and D. Lorenz: Complexity and Applications of the Homotopy
%     Principle for Uniformly Constrained Sparse Minimization
%
% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 05-March-2019

% add paths to data and auxiliary functions
addpath data/mnist_database/;
addpath ..;
warning off;

% read MNIST data
mnistImagesTrain = loadMNISTImages('train-images-idx3-ubyte');
mnistLabelsTrain = loadMNISTLabels('train-labels-idx1-ubyte');
mnistImagesTest = loadMNISTImages('t10k-images-idx3-ubyte');
mnistLabelsTest = loadMNISTLabels('t10k-labels-idx1-ubyte');

% choose two digits for binary classification
digitX = 0;
digitY = 1;

% extract training data
XTrain = mnistImagesTrain(:, mnistLabelsTrain == digitX);
YTrain = mnistImagesTrain(:, mnistLabelsTrain == digitY);

% extract test data
XTest = mnistImagesTest(:, mnistLabelsTest == digitX);
YTest = mnistImagesTest(:, mnistLabelsTest == digitY);

% choose sample sizes
nX = size(XTrain, 2);
nY = size(YTrain, 2);
n = nX + nY;

% calculate homotopy path
XBar = mean(XTrain, 2);
YBar = mean(YTrain, 2);
SigmaHat = ((XTrain - XBar) * (XTrain - XBar)' + (YTrain - YBar) * (YTrain - YBar)') / n;
[betaPathHoudini, lambdaPathHoudini] = l1HoudiniPath(SigmaHat, XBar - YBar, 0);
betaMin = betaPathHoudini(:, end);

% apply path validation algorithm
betaVal = sldaValFunction(XTrain, YTrain, XTest, YTest);

% save images for paper
% XBar
f=figure('Visible','off');
imagesc(reshape(XBar, 28, 28));
axis equal;
set(gca, 'Visible', 'off');
colormap gray;
colorbar('Location', 'manual', 'Units', 'pixels', 'Position', [470 48 25 342]);
saveas(f,'exampleMNISTXBar.png');
close(f);
% YBar
f=figure('Visible','off');
imagesc(reshape(YBar, 28, 28));
axis equal;
set(gca, 'Visible', 'off');
colormap gray;
colorbar('Location', 'manual', 'Units', 'pixels', 'Position', [470 48 25 342]);
saveas(f,'exampleMNISTYBar.png');
close(f);
% deltaHat
f=figure('Visible','off');
imagesc(reshape(XBar - YBar, 28, 28));
axis equal;
set(gca, 'Visible', 'off');
colormap gray;
colorbar('Location', 'manual', 'Units', 'pixels', 'Position', [470 48 25 342]);
saveas(f,'exampleMNISTDeltaHat.png');
close(f);
% SigmaHat
f=figure('Visible','off');
imagesc(reshape(betaMin, 28, 28));
axis equal;
set(gca, 'Visible', 'off');
colormap gray;
colorbar('Location', 'manual', 'Units', 'pixels', 'Position', [470 48 25 342]);
saveas(f,'exampleMNISTBetaMin.png');
close(f);
% betaHat
f=figure('Visible','off');
imagesc(reshape(betaVal, 28, 28));
axis equal;
set(gca, 'Visible', 'off');
colormap gray;
colorbar('Location', 'manual', 'Units', 'pixels', 'Position', [470 48 25 342]);
saveas(f,'exampleMNISTBetaVal.png');
close(f);
% SigmaHat * betaHat
f=figure('Visible','off');
imagesc(reshape(SigmaHat * betaVal, 28, 28));
axis equal;
set(gca, 'Visible', 'off');
colormap gray;
colorbar('Location', 'manual', 'Units', 'pixels', 'Position', [470 48 25 342]);
saveas(f,'exampleMNISTSigmaBetaVal.png');
close(f);