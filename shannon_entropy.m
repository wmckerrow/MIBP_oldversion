function entropy = shannon_entropy(x)
% This function computes the entropy of one or more discrete random variables. It is
% used to find the entropy of all base pairs.
% INPUT: x, a k x m matrix, where k is the number of values the variables
%          can on (so k = 2 in the case of base pairs), and m is the number
%          of variables. Each entry in a column is the probability that the
%          variable takes on the corresponding value.
% OUTPUT: entropy, an m-length array, where the i'th entry is the entropy
%           of the i'th random variable. 

x(x==0) = 1; % a hack to get around NaN
entropy = -sum(x .* log2(x));
end

