%EXAMPLEMNISTONEVSREST - reproduces Fig. 12-13 in (*)
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,
%       dualActiveSet.m, sldaValFunction.m, sldaValHeuristic.m,
%       sldaThreshold.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET,
%           SLDAVALFUNCTION,  SLDAVALHEURISTIC, SLDATHRESHOLD
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

% disable warnings
warning off;

% read MNIST data
mnistImagesTrain = loadMNISTImages('train-images-idx3-ubyte');
mnistLabelsTrain = loadMNISTLabels('train-labels-idx1-ubyte');
mnistImagesTest = loadMNISTImages('t10k-images-idx3-ubyte');
mnistLabelsTest = loadMNISTLabels('t10k-labels-idx1-ubyte');

% select digit for one-vs-rest classification
digit = 5;

% extract training data
XTrain = mnistImagesTrain(:, mnistLabelsTrain == digit);
YTrain = mnistImagesTrain(:, mnistLabelsTrain ~= digit);

% extract test data
XTest = mnistImagesTest(:, mnistLabelsTest == digit);
YTest = mnistImagesTest(:, mnistLabelsTest ~= digit);

% calculate accuracy paths
[~, lambdaVal, ~, plotData] = sldaValFunction(XTrain, YTrain, XTest, YTest);
[~, lambdaHeu, ~, ~, lambdaPathHeu, accPathHeu] = sldaValHeuristic(XTrain, YTrain, XTest, YTest);

% generate plot
close all;
set(gcf, 'visible', 'off');
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [19 5])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition', [0 0 19 5])
set(gca, 'FontSize', 8)
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'Box', 'on')
xlabel('$\lambda, \ \lambda^k$', 'Interpreter', 'latex');
set(gca, 'YLim', [0, 1]);
set(gca, 'XLim', [0, plotData{4}]);
set(gca, 'XLim', [0, .38]);
hold on;
h2 = plot(lambdaPathHeu, accPathHeu, '-g.', 'LineWidth', 1);
h1 = stairs(plotData{1}, plotData{2}, 'Color', 'k', 'LineWidth', 1);
line([lambdaVal, lambdaVal], get(gca, 'YLim'), 'Color', 'k', 'LineStyle', ':')
line([lambdaHeu, lambdaHeu], get(gca, 'YLim'), 'Color', 'g', 'LineStyle', ':')
legend([h1, h2], {'$\varepsilon(\lambda)$', '$\varepsilon(\beta^k, \tau(\beta^k))$'}, 'Interpreter', 'latex', 'Location', 'southwest');
system('sudo rm exampleMNISTOneVsRest.pdf');
print(gcf, 'exampleMNISTOneVsRest', '-dpdf', '-r0');
system('sudo pdfcrop exampleMNISTOneVsRest.pdf exampleMNISTOneVsRest.pdf');
system('evince exampleMNISTOneVsRest.pdf &');