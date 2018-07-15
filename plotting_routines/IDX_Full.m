function [idx,idxOdd,idxEven] = IDX_Full(M)

counter = 1;
for ii = 0 : M
    for jj = 0 : ii
        idx(counter,:) = [M - ii, (ii - jj), jj];
        counter = counter + 1;
    end
end

% reduce to 2D
idx = idx(rem(idx(:,3),2) == 0,:);

% odd variables with respect to x_1
idxOdd = idx(rem(idx(:,1),2) ~= 0,:);

% even variables with respect to x_1
idxEven = idx(rem(idx(:,1),2) == 0,:);

end