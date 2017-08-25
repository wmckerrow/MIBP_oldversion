function [tree_strings,all_tree_strings] = full_split(split_cluster,tree_strings,all_tree_strings,split_bp)
% Split one of the regions into one with a certain base pair and one without it.
% INPUTS: The cluster to split, the names of all the clusters, the names of all the tree nodes, the base pair that will define the split.
% OUTPUTS: The new clusters, the names of the nodes in the new tree.

program_constants;

disp(['splitting ' num2str(split_bp(1)) ' ' num2str(split_bp(2)) ' in cluster ' split_cluster])

split_node_index =  ismember(tree_strings,split_cluster);

%create new mibp files for children of split_cluster
str1 = strcat(split_cluster, '0');
str2 = strcat(split_cluster, '1');
mibp_file = strcat(RNA_NAME, '_', split_cluster, '_mibps.txt');
mibp_file1 = strcat(RNA_NAME, '_', str1, '_mibps.txt');
mibp_file2 = strcat(RNA_NAME, '_', str2, '_mibps.txt');
if strcmp('', split_cluster)
    old_mibps = [];
else
    old_mibps = dlmread(mibp_file);
end
new_mibps = [old_mibps; split_bp];
dlmwrite(mibp_file1, new_mibps);
dlmwrite(mibp_file2, new_mibps);

%split the cluster
B_split_cluster(split_cluster, old_mibps, split_bp);

%create and save MI of two new clusters
% [max_mil max_pairl] = B_get_MI(str1);
% [max_mir max_pairr] = B_get_MI(str2);
% left_MI = strcat(RNA_NAME, '_', str1, '_MI.txt');
% right_MI = strcat(RNA_NAME, '_', str2, '_MI.txt');
% dlmwrite(left_MI, [max_mil max_pairl]);
% dlmwrite(right_MI, [max_mir max_pairr]);

tree_strings(split_node_index) = []; % eliminate the split cluster from tree_strings
all_tree_strings{end+1} = str1;
all_tree_strings{end+1} = str2;
tree_strings{end+1} = str1;
tree_strings{end+1} = str2;


make_all_probs_files_vB({str1,str2});

end