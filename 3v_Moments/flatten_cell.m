% given an array of cells, the following routine converts it to a vector
function f = flatten_cell(data)
f = [];

for i = 1 : length(data)
    for j = 1 : length(data{i})
        f = [f data{i}(j)];
    end
end
end