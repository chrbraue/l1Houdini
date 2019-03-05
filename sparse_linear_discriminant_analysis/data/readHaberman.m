% read data
data = csvread('haberman.dat');
features = data(:, 1:end-1);
labels = data(:, end);

% split data into two classes
X = features(labels == 1, :)';
Y = features(labels == 2, :)';

% split data into training and validation sets
nX = size(X, 2);
nY = size(Y, 2);
pX = randperm(nX);
pY = randperm(nY);
XTrain = X(:, pX(1:ceil(nX / 2)));
YTrain = Y(:, pY(1:ceil(nY / 2)));
XVal = X(:, pX(ceil(nX / 2) + 1:end));
YVal = Y(:, pY(ceil(nY / 2) + 1:end));