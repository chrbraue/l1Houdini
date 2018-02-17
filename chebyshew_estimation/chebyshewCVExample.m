%CHEBYSHEWCVEXAMPLE - random example for chebyshewCVFunction
%
% This example reproduces Fig. 2 in [Brauer et al. 2018]
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,
%       dualActiveSet.m, chebyshewCVFunction.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET,
% CHEBYSHEWCVFUNCTION

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 17-February-2018

% add path to algorithm
addpath ..;

% set seed value
rng(1);

% choose m, n, ground truth sparsity and magnitude of noise
foldSize = 10;
N = 10;
m = N * foldSize;
n = 50;
sparsity = 5;
deltaGT = 1;

% generate A, ground truth (normally distributed) and noise (uniformly
% distributed)
A = randn(m, n);
xGT = randn(n, 1);
p = randperm(n);
xGT(p(sparsity + 1:end)) = 0;
eta = deltaGT * (2 * rand(m, 1) - 1);

% calculate measurement vector b
b = A * xGT + eta;

% perform cross-validation and compute distance of xCV to ground truth
[xCV, deltaCV, errCV, pathDeltaCV, pathErrorCV] = chebyshewCVFunction(A, b, N);
distCV = norm(xCV - xGT, 1);

% solve with known/true delta and compute distance to ground truth
xDeltaGT = l1Houdini(A, b, deltaGT);
distGT = norm(xDeltaGT - xGT, 1);

% compute full solution path and find the particular kink with smallest
% distance to ground truth
[pathXFull, pathDeltaFull] = l1HoudiniPath(A, b, 0);
pathXFullDist = sum(abs(pathXFull - xGT));
[distBest, index] = min(pathXFullDist);
deltaBest = pathDeltaFull(index);
xBest = pathXFull(:, index);

% compare estimates w.r.t. distance to ground truth
fprintf('\ndeltaGT = %f\n', deltaGT);
fprintf('||xDeltaGT - xGT||_1 = %f\n', distGT);
fprintf('\ndeltaCV = %f\n', deltaCV);
fprintf('||xCV - xGT||_1 = %f\n', distCV);
fprintf('\ndeltaBest = %f\n', deltaBest);
fprintf('||xBest - xGT||_1 = %f\n', distBest);
fprintf('\n');

% interpolate error and distance
pathUnion = union(pathDeltaCV, pathDeltaFull);
pathUnion = sort(pathUnion, 'descend');
pathErrorMean = interp1(pathDeltaCV, pathErrorCV, pathUnion);
pathXFullDist = interp1(pathDeltaFull, pathXFullDist, pathUnion);

% initialize figure
figure;
set(gca, 'TickLabelInterpreter', 'latex');
hold on;

% plot average testing error and distance to ground truth
h1 = semilogy(pathUnion(1) - pathUnion, pathErrorMean, 'k--');
h2 = semilogy(pathUnion(1) - pathUnion, pathXFullDist, 'k-');

% set/initialize axis/tick labels
xlabel('$\delta$', 'Interpreter', 'latex');
yTicks = [];
yTickLabels = {};
xTicks = [];
xTickLabels = {};

% set horizontal axis limit
set(gca, 'XLim', [0, pathUnion(1) - pathUnion(end)]);

% plot horizontal line at minimum cv error
[~, index] = min(pathErrorMean);
line(get(gca, 'XLim'), [pathErrorMean(index) pathErrorMean(index)], 'Color', [.5 .5 .5], 'LineWidth', .1);
yTicks = [yTicks, pathErrorMean(index)];
yTickLabels{end + 1} = sprintf('$\\varepsilon_{\\mathrm{cv}}$');

% plot horizontal line at distance of xCV to ground truth
line(get(gca, 'XLim'), [distCV distCV], 'Color', [.5 .5 .5], 'LineWidth', .1);
yTicks = [yTicks, distCV];
yTickLabels{end + 1} = sprintf('$d_{\\mathrm{cv}}$');

% plot horizontal line at distance of xDeltaGT to ground truth
line(get(gca, 'XLim'), [distGT distGT], 'Color', [.5 .5 .5], 'LineWidth', .1, 'LineStyle', '--');
yTicks = [yTicks, distGT];
yTickLabels{end + 1} = sprintf('$d^{\\dagger}$');

% adapt limits of vertical axis
yLim = get(gca, 'YLim');
set(gca, 'YLim', [yLim(1) / 1.1, yLim(2) * 1.1]);

% set vertical ticks and tick labels
yTicks = [yTicks, max(pathXFullDist)];
yTickLabels{end + 1} = sprintf('$%2.2f$', max(pathXFullDist));
[yTicks, p] = sort(yTicks);
yTickLabels = yTickLabels(p);
set(gca, 'YTick', yTicks, 'YTickLabel', yTickLabels);

% plot vertical line at cv delta
line(pathDeltaCV(1) - [deltaCV deltaCV], get(gca, 'YLim'), 'Color', [.5 .5 .5], 'LineWidth', .1);
xTicks = [xTicks pathDeltaCV(1) - deltaCV];
xTickLabels{end + 1} = sprintf('$\\delta_{\\mathrm{cv}}$');

% plot vertical line at original delta
line(pathDeltaCV(1) - [deltaGT deltaGT], get(gca, 'YLim'), 'Color', [.5 .5 .5], 'LineWidth', .1, 'LineStyle', '--');
xTicks = [xTicks pathDeltaCV(1) - deltaGT];
xTickLabels{end + 1} = sprintf('$\\delta^{\\dagger}$');

% set horizontal ticks and tick labels
xTicks = [xTicks, 0, pathUnion(1) - pathUnion(end)];
xTickLabels{end + 1} = sprintf('$%2.2f$', max(pathUnion));
xTickLabels{end + 1} = sprintf('$%2.2f$', min(pathUnion));
[xTicks, p] = sort(xTicks);
xTickLabels = xTickLabels(p);
set(gca, 'XTick', xTicks, 'xTickLabel', xTickLabels);

% set legend
legend([h1, h2], {'$\varepsilon(\delta)$', '$||x^*(\delta) - x^{\dagger}||_1$'}, 'Interpreter', 'latex');