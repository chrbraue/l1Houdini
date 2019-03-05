function A = readSparseInteger(fileName, p, n)

fid = fopen(fileName);

indexValuePairs = textscan(fid, '%d:%d');

i = 1;

A = zeros(p, n);

prevIndex = 1;

for k = 1:length(indexValuePairs{1})    
    index = indexValuePairs{1}(k);
    value = indexValuePairs{2}(k);
    if index < prevIndex
        i = i + 1;
    end
    prevIndex = index;
    A(index, i) = value;
end

end