function [forced_pairs,forced_nonpairs] = B_entropy_constraints3(tree_path, entropy_cutoff2, pairs)
% This function applies the "entropy constraints." That is, it creates a
% new secondary structure probability distribution where base pairs with
% entropy below "entropy_cutoff" are forced to be paired or unpaired,
% depending on whether the pairing probability is close to 0 or 1.
% INPUTS: tree_path, a string of 0's and 1's identifying the cluster
%         entropy_cutoff, a double between 0 and 1, the entropy at which a
%           base pair is forced to pair or not pair
% OUTPUT: Q, a double, the partition function of the constrained model.
%           (This function also write several files.)

program_constants;

% get bp probability file. it's an mx3 matrix, where each row has form:
% [base 1, base 2, probability of pairing]
%pfsfile = strcat('B_', RNA_NAME, '_', tree_path, '.pfs');
probs_file = strcat(RNA_NAME, '_', tree_path, '_probs.txt');
%system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PROB_PROG ' -t ' pfsfile ' ' probs_file ' >/dev/null']);
bp_probs = read_bp_probs(probs_file);

bp_probs = bp_probs(bp_probs(:,3) > 0, :); % eliminate pairs with probability 0
bp_entropies = bp_probs; % we will change column 3 of bp_entropies from probabilities to entropies
prob_pairs = [bp_probs(:,3), 1-bp_probs(:,3)]';
entropies = shannon_entropy(prob_pairs)';
for i = 1:size(bp_entropies, 1)
    p = prob_pairs(:,i);
    bp_entropies(i,3) = shannon_entropy(p);
end


% find base pairs with low entropy and constrain them to pair or nonpair:

len_bp_probs = size(bp_probs, 1);
k = 1; % k is the index of the base pair in bp_probs. start with first base pair and iterate through list.
forced_pairs = pairs;
forced_nonpairs = [];
for i = 1:RNA_LENGTH % two loops iterate through all pairs of bases
   for j = (i+1):RNA_LENGTH
      if (k > len_bp_probs || bp_probs(k,1) ~= i || bp_probs(k,2) ~= j) %if current (i,j) pair doesn't equal k'th nuc pair in bp_probs
          % in this case, (i,j) cannot pair, and so
          % probability of pairing is 0 and so is its entropy
          
          forced_nonpairs = [forced_nonpairs; [i,j]];
          
      elseif (bp_entropies(k, 3) < entropy_cutoff2) %if current (i,j) pair is in bp_probs and has low entropy

          if (bp_probs(k, 3) > 0.5) % (i,j) is a forced pair
              if isempty(pairs) || ~ismember([i,j],pairs,'rows')
                forced_pairs = [forced_pairs; [i, j]];
              end
          else % or else it's a forced nonpair
              forced_nonpairs = [forced_nonpairs; [i, j]];
          end
          k = k + 1;
      else % else current (i,j) pair is the k'th pair in bp_probs, and it has entropy >= entropy_cutoff
          k = k + 1;
      end
   end
end

forced_pairs = sortrows(forced_pairs);

end

