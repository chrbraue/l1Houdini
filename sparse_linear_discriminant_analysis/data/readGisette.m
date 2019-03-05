% read data
trainingData = dlmread('data/nips2003_fsc/GISETTE/gisette_train.data');
trainingLabels = dlmread('data/nips2003_fsc/GISETTE/gisette_train.labels');
validationData = dlmread('data/nips2003_fsc/GISETTE/gisette_valid.data');
validationLabels = dlmread('data/nips2003_fsc/GISETTE/gisette_valid.labels');

% split data
XTrain = trainingData(trainingLabels == 1, :)';
YTrain = trainingData(trainingLabels == -1, :)';
XVal = validationData(validationLabels == 1, :)';
YVal= validationData(validationLabels == -1, :)';