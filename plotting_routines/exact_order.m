% x1 is the first x coordinate
% y1 is the first y coordinate
% x2 is the final x coordinate
% order is the desired order needed 
function[reference_order] = exact_order(x1,y1,x2,order)

y2 = y1 * exp(order * log(x2/x1));
reference_order = [x1 x2 y1 y2];
end