function [tree_strings,all_tree_strings] = add_hibp(tree_strings,all_tree_strings)
% Keep imposing high probability base pairs until there are non left.
% INPUTS: Node labels for regions and internal tree nodes.
% OUTPUTS: Updated node labels for regions and internal tree nodes.

program_constants;
old_tree_strings = tree_strings;
for i = 1:length(old_tree_strings) % Loop through all clusters
    cluster = old_tree_strings{i};
    % Need to keep track of forced pairs so that we can avoid pairs that
    % are in a common stem.
    if isempty(cluster)
        forced_pairs = [];
    else
        mibps = dlmread([RNA_NAME '_' cluster '_mibps.txt']);
        forced_pairs = mibps(cluster == '1',:);
    end
    total_bp_list=hibpstems(cluster,forced_pairs);
    %List all high probability base that are not in the same stem as a
    %pair that is already forced. Keep adding until there are none left.
    while ~isempty(total_bp_list)
        disp(['hibp split at ' cluster])
        % Split into a cluster with and without the high probability pair
        [tree_strings,all_tree_strings] = full_split(cluster,tree_strings,all_tree_strings,total_bp_list(1,1:2));
        % Move to the new cluster and look for another high probability
        % pair.
        cluster = [cluster,'1'];
        mibps = dlmread([RNA_NAME '_' cluster '_mibps.txt']);
        forced_pairs = mibps(cluster == '1',:);
        total_bp_list=hibpstems(cluster,forced_pairs);
    end
end

end