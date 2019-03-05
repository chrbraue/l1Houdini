function A = readSparseBinary(fileName, p, n)

fid = fopen(fileName);

indices = textscan(fid, '%d');

i = 1;

A = zeros(p, n);

prevIndex = 1;

for k = 1:length(indices{1})    
    index = indices{1}(k);
    if index < prevIndex
        i = i + 1;
    end
    prevIndex = index;
    A(index, i) = 1;
end

end