% read data
trainingData = dlmread('data/nips2003_fsc/ARCENE/arcene_train.data');
trainingLabels = dlmread('data/nips2003_fsc/ARCENE/arcene_train.labels');
validationData = dlmread('data/nips2003_fsc/ARCENE/arcene_valid.data');
validationLabels = dlmread('data/nips2003_fsc/ARCENE/arcene_valid.labels');

% split data
XTrain = trainingData(trainingLabels == 1, :)';
YTrain = trainingData(trainingLabels == -1, :)';
XVal = validationData(validationLabels == 1, :)';
YVal= validationData(validationLabels == -1, :)';