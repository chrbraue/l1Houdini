% plotMNIST

fprintf('data import...');
X = loadMNISTImages('train-images-idx3-ubyte');
fprintf(' finished\n');
shape = [2, 5];
padding = 0;
m = size(X, 2);
n = sqrt(size(X, 1));
indices= [2 4 6 8 3 1 14 16 18 5];
clear tmp;
k = 1;
for i = 0:shape(1) - 1    
    for j = 0:shape(2) - 1
        tmp(i * (n + padding) + 1:i * (n + padding) + n, ...
            j * (n + padding) + 1:j * (n + padding) + n) ...
            = reshape(1 - X(:, indices(k)), 28, 28);
        k = k + 1; 
    end
end
tmp = [zeros(padding, size(tmp, 2)); tmp; zeros(padding, size(tmp, 2))];
tmp = [zeros(size(tmp, 1), padding), tmp, zeros(size(tmp, 1), padding)];
imwrite(tmp, 'mnistDigits.png');
