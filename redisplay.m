pg_eval; % set program constants, such as RNA_NAME, FORCED_PAIRS, FORCED_NONPAIRS, etc.

%% Save everything
load([RNA_NAME '.mat']);

disp(['native_clust_prob:',num2str(probs(find(ismember(tree_strings,native_cluster))))])
disp(['exp_clust_prob:',num2str(sum(probs.^2))])
disp(['rand_clust_prob:',num2str(mean(probs))])
disp(['num_structs:',num2str(sum(structs)),',',num2str(sum(structs_le))])
disp(['prob_le:',num2str(sum(probs_le))])
disp(['entropy:',num2str(sum_of_bp_entropy),',',num2str(sum(probs.*entropies))])


%% Clean up
system('rm temp*');
