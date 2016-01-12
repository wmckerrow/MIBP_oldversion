function total_bp_list=hibpstems(clus,forced_pairs)
% Find probable basepairs not already in the same stem as one of the forced
% pairs.
% INPUTS: The name of the region of interest. A list of all the pairs
% forced in that region.
% OUTPUTS: A list of all highly probable base pairs not in a stem that
% already has a forced pair.

total_bp_list = [];

program_constants;

% Read the probability of all the base pairs.
prob_file = strcat(RNA_NAME, '_',clus,'_probs.txt');
base_pairs = read_bp_probs(prob_file);

% Find all bps in the same stem as a forced pair
forced_stems = [0 0];
for j = 1:size(forced_pairs,1)
    forced_stems=[forced_stems;find_stem_bp(forced_pairs(j,:),clus,stem_thresh)];
end

% Possible high probability bps must have probability abouve the
% hibp_cutoff
eligible_list = base_pairs(base_pairs(:,3) > hibp_cutoff,:);

% Check whether an eligible base pair is in the forced stems list. If it is
% not add it
add_here = 1;
for j = 1:size(eligible_list,1)
    [~,next_largest] = max(eligible_list(:,3));
    pair = eligible_list(next_largest,:);
    eligible_list(next_largest,:) = [];
    if ~ismember(pair(1:2),forced_stems,'rows')
        total_bp_list(add_here,:) = pair;
        forced_stems=[forced_stems;find_stem_bp(pair(1:2),clus,stem_thresh)];
        add_here = add_here+1;
    end
end

end