% read data
trainingData = readSparseInteger('data/nips2003_fsc/DEXTER/dexter_train.data', 20000, 300);
trainingLabels = dlmread('data/nips2003_fsc/DEXTER/dexter_train.labels');
validationData = readSparseInteger('data/nips2003_fsc/DEXTER/dexter_valid.data', 20000, 300);
validationLabels = dlmread('data/nips2003_fsc/DEXTER/dexter_valid.labels');

% split data
XTrain = trainingData(:, trainingLabels == 1);
YTrain = trainingData(:, trainingLabels == -1);
XVal = validationData(:, validationLabels == 1);
YVal= validationData(:, validationLabels == -1);