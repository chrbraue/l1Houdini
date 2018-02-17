function [xCV, deltaCV, errCV, pathDelta, pathError] = chebyshewCVFunction(A, b, N)
%CHEBYSHEWCVFUNCTION - cross validation scheme for chebyshew estimation
%chebyshewCVFunction(A, b, N) performs an N-fold cross-validation scheme
%for chebyshew estimation using the l1-Houdini algorithm
%
% Syntax: [xCV, deltaCV, errCV, pathDelta, pathError] = chebyshewCVFunction(A, b, N)
%
% Input:
%    A - 
%    b - 
%    N - number of folds (integer)
%    
% Output:
%    xCV - linear predictor (n vector)
%    deltaCV - tuning parameter associated wit xCV
%    errCV - predicted error
%    pathDelta - all parameters associated with breakpoints
%    pathError - all errors associated with these breakpoints
%
% Example: see chebyshewCVExample.m
%
% Other m-files required: l1HoudiniPath.m, primalActiveSet.m,  dualActiveSet.m
% Subfunctions: none
% MAT-files required: none
%
% See also: L1HOUDINIPATH, PRIMALACTIVESET,  DUALACTIVESET

% Author: Christoph Brauer (TU Braunschweig)
% contact email address: ch.brauer@tu-braunschweig.de
% February 2018; Last revision: 17-February-2018

% initialize variables
pathX = {};
pathsDelta = {};
pathsError = {};
pathDelta = [];

% get number of samples
m = size(A, 1);

% generate random folds via permuation
p = randperm(m);

% compute fold size
foldSize = floor(m / N);

% solve subproblems
for k = 1:N
    % determine testing and training folds
    testingFold = p((k - 1) * foldSize + 1:(k - 1) * foldSize + foldSize);
    trainingFold = setdiff(1:m, testingFold);
    % solve k-th subproblem
    [pathX{k}, pathsDelta{k}] = l1HoudiniPath(A(trainingFold, :), b(trainingFold), 0); %#ok<*AGROW>
    % compute testing errors along the entire path
    pathsError{k} = max(abs(A(testingFold, :) * pathX{k} - b(testingFold)));
    % add breakpoints to final parameter path
    pathDelta = union(pathDelta, pathsDelta{k});
end
% sort final parameter path in descending order
pathDelta = sort(pathDelta, 'descend');

% interpolate missing error values and calculate mean error
pathError = zeros(size(pathDelta));
for k = 1:N
    interp = interp1(pathsDelta{k}, pathsError{k}, pathDelta);
    % the k-th solution path is constant for delta >= max(pathDelta{k}) and
    % the corresponding error is pathError{k}(1)
    interp(isnan(interp) & pathDelta >= max(pathsDelta{k})) = pathsError{k}(1);
    % the k-th solution is not defined for delta < delta^k_min, so we
    % remove the corresponding entries from the final parameter path
    index = ~isnan(interp);
    interp = interp(index);
    pathDelta = pathDelta(index);
    % update mean error
    pathError = pathError(index) + interp / N;
end

% solve again with complete data and optimal parameter
[errCV, index] = min(pathError);
deltaCV = pathDelta(index);
xCV = l1Houdini(A, b, deltaCV);