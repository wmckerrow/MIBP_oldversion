function [max_mi max_pair sum_mi] = rna_mi3(samples)
% This function finds the mutual information (MI) of all pairs of nucleotides
% from a sampled collection of secondary structures.
% INPUT: structures, a matrix of base pairs (see readct.m)
% OUTPUTS: max_mi, a num, the maximum average MI of any pair of nucleotides
%            with all other pairs of nucleotides.
%          max_pair, a length-2 array of ints, the pair of nucleotides
%            achieving the max average MI
%          avg_mi, an m x 3 matrix, where the first two entries of a row
%            are a pair of nucleotides, and the third entry is the average
%              MI this pair has with all other pairs.

eta = 0.01;  % throw away pairs with prob less than eta or greater than 1-eta
             % when calculating MIBP (see luan lin's thesis for why) 
delta = 0.001;  % minimum bp prob

counts = count_from_struct(samples); % make a table of base pairs for each structure

base_pairs = pair_probs2(counts);  %  get the base pair probs
base_pairs = base_pairs(base_pairs(:,3) > delta,:); % filter out pairs with low prob
% (see Luan Lin's thesis for theoretical justification of filtering out low
% probs)

% do the real work here
% loop through all pairs of base pairs
% this can probably be made faster
mi = zeros(length(base_pairs), length(base_pairs));
for i=1:(size(base_pairs,1)-1)
    if base_pairs(i,3) > eta && base_pairs(i,3) < 1-eta   % get rid of low entropy pairs
        %calculate the mutual information between base pair i and base pair j
        p1_count = samples(:,base_pairs(i,1)) == base_pairs(i,2); % get a vector of all structures that have this bp
        for j=(i+1):size(base_pairs,1)            
            if base_pairs(j,3) > eta && base_pairs(j,3) < 1-eta
                p2_count = samples(:,base_pairs(j,1)) == base_pairs(j,2);
                
                mi(i,j) = mutual_info(p1_count, p2_count); % get the MI between these two bps using the structures
                mi(j,i) = mi(i,j);  % MI is symmetric
            end
        end;
    end
end

% get the max average MI, throwing away pairs with low entropy
% for why we throw away pairs with low entropy, see Luan Lin's thesis
c = 1;
sum_mi = zeros(size(base_pairs,1), 3);
for i=1:length(base_pairs)
    if base_pairs(i,3) < 1-eta && base_pairs(i,3) > eta
        sum_mi(c, 1:2) = base_pairs(i, 1:2);
        sum_mi(c,3) = sum(mi(i,:));
        c = c + 1;
    end
end
sum_mi(sum_mi(:,3) == 0, :) = [];

[max_mi ndx] = max(sum_mi(:,3));
max_pair = sum_mi(ndx, 1:2);

end

