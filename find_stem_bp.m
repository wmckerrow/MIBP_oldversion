function stem = find_stem_bp(pair,cluster,stem_variation)
% Finds other base pairs in the same 'stem' by following adjacent pairs in the forward direction and backward direction until the probability falls below a threshold.
% INPUTS: The pair that we want to find the stem for, the region (conditional distribution) that we are considering, the amount we allow probabilities to vary in a single stem (usually 2 fold)

program_constants;

% Read base pair probabilities and put into into an easy to access array.
cluster_prob_file = strcat(RNA_NAME, '_', cluster, '_probs.txt');
cluster_probs = read_bp_probs(cluster_prob_file);
cluster_prob_mat = zeros(RNA_LENGTH,RNA_LENGTH);
for i = 1:length(cluster_probs)
    cluster_prob_mat(cluster_probs(i,1),cluster_probs(i,2))=cluster_probs(i,3);
    cluster_prob_mat(cluster_probs(i,2),cluster_probs(i,1))=cluster_probs(i,3);
end

% Look for adjacent base pairs until probability falls below this threshold.
threshold = cluster_prob_mat(pair(1),pair(2))/stem_variation;

stem=pair;

% Look for neighbors in the forward direction
next1 = mod(pair(1),RNA_LENGTH)+1;
next2 = mod(pair(2)-2,RNA_LENGTH)+1;
next_prob = cluster_prob_mat(next1,next2);
while next_prob > threshold % Stop when we fall below the threshold
    stem = [stem;next1,next2];
    next1 = mod(next1,RNA_LENGTH)+1;
    next2 = mod(next2-2,RNA_LENGTH)+1;
    next_prob = cluster_prob_mat(next1,next2);
end

%Look in the backward direction
prev1 = mod(pair(1)-2,RNA_LENGTH)+1;
prev2 = mod(pair(2),RNA_LENGTH)+1;
prev_prob = cluster_prob_mat(prev1,prev2);
while prev_prob > threshold
    stem = [stem;prev1,prev2];
    prev1 = mod(prev1-2,RNA_LENGTH)+1;
    prev2 = mod(prev2,RNA_LENGTH)+1;
    prev_prob = cluster_prob_mat(prev1,prev2);
end

end