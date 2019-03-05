% read data
data = csvread('processed.cleveland-clean.dat');
features = data(:, 1:end-1);
labels = data(:, end);

% split data into two classes
X = features(labels == 0, :)';
Y = features(labels >= 1, :)';

% split data into training and validation sets
nX = size(X, 2);
nY = size(Y, 2);
pX = randperm(nX);
pY = randperm(nY);
XTrain = X(:, pX(1:80));
YTrain = Y(:, pY(1:80));
XVal = X(:, pX(80 + 1:end));
YVal = Y(:, pY(80 + 1:end));