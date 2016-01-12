function [ base_pairs ] = pair_probs2(counts)
% Converts the counts returned by count_pairs.m into base pair
% probabilities.
% INPUT: structure, an L x L array (see readct.m) recording base
%          pairs
% OUTPUT: base_pairs, an m x 3 matrix, where m is the number of distinct base pairs
%           that were observed in the sampled structures. base_pairs(j, 1)
%           and base_pairs(j, 2) are the nucleotides of the j'th pair of
%           bases. nucleotides(j, 3) is the sample probability of this base
%           pair forming

program_constants;

[r c] = find(counts>0);
base_pairs = zeros(length(r),3);
base_pairs(:,1) = r;
base_pairs(:,2) = c;
base_pairs(:,3) = counts(counts>0)/NUM_SAMPLES;

base_pairs = sortrows(base_pairs, 1:2);
end

