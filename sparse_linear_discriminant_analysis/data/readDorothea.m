% read data
trainingData = readSparseBinary('data/nips2003_fsc/DOROTHEA/dorothea_train.data', 100000, 800);
trainingLabels = dlmread('data/nips2003_fsc/DOROTHEA/dorothea_train.labels');
validationData = readSparseBinary('data/nips2003_fsc/DOROTHEA/dorothea_valid.data', 100000, 800);
validationLabels = dlmread('data/nips2003_fsc/DOROTHEA/dorothea_valid.labels');

% split data
XTrain = trainingData(:, trainingLabels == 1);
YTrain = trainingData(:, trainingLabels == -1);
XVal = validationData(:, validationLabels == 1);
YVal= validationData(:, validationLabels == -1);