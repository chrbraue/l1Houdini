%EXAMPLEMNISTONEVSONE - reproduces Fig. 4-5 in (*)
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

% choose plot colors
plotColors = {'b', 'g', 'k'};

% loop over all digit pairs
digitX = 4;
digitY = 9;

% extract training data
XTrain = mnistImagesTrain(:, mnistLabelsTrain == digitX);
YTrain = mnistImagesTrain(:, mnistLabelsTrain == digitY);

% extract test data
XTest = mnistImagesTest(:, mnistLabelsTest == digitX);
YTest = mnistImagesTest(:, mnistLabelsTest == digitY);

% choose sample sizes
nX = [5, 25, size(XTrain, 2)];
nY = [5, 25, size(YTrain, 2)];
n = nX + nY;

% initialize figure
close all;
set(gcf, 'visible', 'off');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [19 7.5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 19 7.5]);
set(gca, 'FontSize', 8);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'Box', 'on');
xlabel('$\lambda$', 'Interpreter', 'latex');
ylabel('$\varepsilon(\lambda)$', 'Interpreter', 'latex');
set(gca, 'YLim', [.3, 1]);
lambdaMax = 1;
set(gca, 'XLim', [0, lambdaMax]);
hold on;
legendInfo = {};
plots = [];

% run path validation algorithm multiple times and plot accuracies
for k = 1:3
    [beta, lambda, accuracy, plotData] = sldaValFunction(XTrain(:, 1:nX(k)), YTrain(:, 1:nY(k)), XTest, YTest);   
    plotData{1}(end) = lambdaMax;                        
    stairs(plotData{1}, plotData{2}, 'Color', plotColors{k}, 'LineWidth', 1);
    plots = [plots stairs(plotData{1}, plotData{2}, 'Color', plotColors{k}, 'LineWidth', 1)];  %#ok<AGROW>
    line([lambda, lambda], get(gca, 'YLim'), 'Color', plotColors{k}, 'LineStyle', ':');
    legendInfo{k} = sprintf('$n = %d$', n(k)); %#ok<SAGROW>
end
legend(plots, legendInfo, 'Interpreter', 'latex', 'Location', 'southwest');
        
% save figure
system('sudo rm exampleMNISTOneVsOne.pdf');
print(gcf, 'exampleMNISTOneVsOne', '-dpdf', '-r0');
system('sudo pdfcrop exampleMNISTOneVsOne.pdf exampleMNISTOneVsOne.pdf');
system('evince exampleMNISTOneVsOne.pdf &');          