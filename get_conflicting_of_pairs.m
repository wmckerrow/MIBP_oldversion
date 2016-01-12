function [conflicting_bp,prob] = get_conflicting_of_pairs(pairs,cluster)
% For a given set base pairs and a region find the base pairs in that region that conflict with all the given pairs and has highest probability. (i.e. the most informative conflicting.)
% INPUTS: the given pair, the region of interest.
% OUTPUTS: base1, base2, probability of pairing for the most informative conflicting.

program_constants;

% Read the pairing proabilities from file.
prob_name = strcat(RNA_NAME, '_',cluster,'_probs.txt');
bp_probs = read_bp_probs(prob_name);

% Make a logical array that is true for conflicting pairs.
conflicting = (bp_probs(:,1) < pairs(1,1) & bp_probs(:,2) > pairs(1,1) & bp_probs(:,2) < pairs(1,2)) | (bp_probs(:,1) > pairs(1,1) & bp_probs(:,1) < pairs(1,2) & bp_probs(:,2) > pairs(1,2)) | bp_probs(:,1) == pairs(1,1) | bp_probs(:,1) == pairs(1,2) | bp_probs(:,2) == pairs(1,1) | bp_probs(:,2) == pairs(1,2);
for i = 2:size(pairs,1)
    conflicting = conflicting & ((bp_probs(:,1) < pairs(i,1) & bp_probs(:,2) > pairs(i,1) & bp_probs(:,2) < pairs(i,2)) | (bp_probs(:,1) > pairs(i,1) & bp_probs(:,1) < pairs(i,2) & bp_probs(:,2) > pairs(i,2)) | bp_probs(:,1) == pairs(i,1) | bp_probs(:,1) == pairs(i,2) | bp_probs(:,2) == pairs(i,1) | bp_probs(:,2) == pairs(i,2));
end

% Among the conflicting base pairs find the one that is most probable.
conf_bp_probs = bp_probs(conflicting,:);
[prob,I] = max(conf_bp_probs(:,3));
conflicting_bp = conf_bp_probs(I,1:2)

end
