%EXAMPLEFSC - reproduces Fig. 7-11 in (*)
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,
%       dualActiveSet.m, sldaValFunction.m, sldaCVFunction.m,
%       sldaValHeuristic.m, sldaThreshold.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET,
%           SLDAVALFUNCTION, SLDACVFUNCTION, SLDAVALHEURISTIC,
%           SLDATHRESHOLD
%
% (*) C. Brauer and D. Lorenz: Complexity and Applications of the Homotopy
%     Principle for Uniformly Constrained Sparse Minimization
%
% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 05-March-2019

% add paths to l1Houdini, data and auxiliary functions
addpath data/;
addpath ../;
warning off;

% choose one of readArcene/readDexter/readDorothea/readGisette/readMadelon/
% readBanknote/readBreast/readCleveland/readHaberman
readBanknote;
% all data sets are from the NIPS 2003 feature selection challenge and
% publicly available under https://archive.ics.uci.edu/ml/index.php

% set string for filename
fileStr = 'exampleFSC';

% set seed value (partition into folds in random)
rng(666, 'twister');

% perform cross validation on training set
N = 10;
[~, lambdaCV, accCV, plotDataCV] = sldaCVFunction(XTrain, YTrain, N);
plotDataCV{5} = [plotDataCV{5}; plotDataCV{5}(end, :)];

% perform path validation using validation set
[~, lambdaVal, accVal, plotDataVal] = sldaValFunction(XTrain, YTrain, XVal, YVal);
[~, lambdaHeu, accHeu, ~, lambdaPathHeu, accPathHeu] = sldaValHeuristic(XTrain, YTrain, XVal, YVal);

% plot results
close all;
set(gcf, 'visible', 'off');
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [19 7.5])
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition', [0 0 19 7.5])
set(gca, 'FontSize', 8)
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'Box', 'on')
xlabel('$\lambda$', 'Interpreter', 'latex');
lambdaMin = plotDataVal{3};
lambdaMax = plotDataCV{4};
% set(gca, 'XLim', [70, 650]); % was used for the Gisette Data Set
set(gca, 'XLim', [lambdaMin, lambdaMax + (lambdaMax - lambdaMin) / 10]);
set(gca, 'YLim', [min(min(plotDataCV{5})) - .025, 1.025]);
hold on;
for k = 1:N
    cvLPlot = stairs(plotDataCV{1}, plotDataCV{5}(:, k), 'Color', [190 190 190] / 256, 'LineWidth', .25);
end
valPlot = stairs([plotDataVal{1}; plotDataCV{1}(end)], [plotDataVal{2}; plotDataVal{2}(end)], 'Color', 'g', 'LineWidth', 1);
line([lambdaVal, lambdaVal], get(gca, 'YLim'), 'Color', 'g', 'LineStyle', ':')
cvPlot = stairs(plotDataCV{1}, plotDataCV{2}, 'k', 'LineWidth', 1);
line([lambdaCV, lambdaCV], get(gca, 'YLim'), 'Color', 'k', 'LineStyle', ':')
heuPlot = plot(lambdaHeu, accHeu, 'kx', 'MarkerSize', 6, 'LineWidth', .5);
legendInfo = {'$\mathrm{CV}(\lambda)$', '$\mathrm{CV}_l(\lambda)$', '$\varepsilon(\lambda)$'};
legend([cvPlot, cvLPlot, valPlot], legendInfo, 'Interpreter', 'latex');

% save figure
system(sprintf('sudo rm %s.pdf', fileStr));
print(gcf, fileStr, '-dpdf', '-r0');
system(sprintf('sudo pdfcrop %s.pdf', fileStr));
system(sprintf('sudo mv %s-crop.pdf %s.pdf', fileStr, fileStr));
system(sprintf('evince %s.pdf &', fileStr));