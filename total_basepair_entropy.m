function total_ent = total_basepair_entropy(tree_path)
% Finds the sum of basepair entropy for the specified cluster.
% sum_X H(X|C=c)
% INPUT: The name of the region of interest.
% OUTPUT: The total base pair entropy in the region

program_constants;

probs_file = strcat(RNA_NAME, '_', tree_path, '_probs.txt');
bp_probs = read_bp_probs(probs_file);
bp_probs = bp_probs(bp_probs(:,3) > 0, :); % eliminate pairs with probability 0
prob_pairs = [bp_probs(:,3), 1-bp_probs(:,3)]';
entropies = shannon_entropy(prob_pairs)';

total_ent = sum(entropies);

end