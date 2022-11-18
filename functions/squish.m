function y = squish(x)
% reshape without argument and remove nans
ss = size(x);
y = reshape(x, prod(ss),1);
y = y(~isnan(y));