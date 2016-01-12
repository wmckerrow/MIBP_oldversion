function B_split_cluster(tree_path, mibps, mibp)
% This function splits a cluster of the tree into its two children. It
% creates the constraint, pfs, ct, and mibp files of the children clusters.
% INPUTS: tree_path, a string of 0's and 1's storing the identity of the
%           cluster to be split
%         mibps, the most informative base pairs in tree_path
%         mibp, the new most informative base pair used to split the
%         cluster

program_constants;

new_mibps = [mibps; mibp];
 

% create left and right tree paths, mibp files, pfs filenames
left_tree_path = strcat(tree_path, '0');
right_tree_path = strcat(tree_path, '1');
left_mibpfile = strcat(RNA_NAME, '_', left_tree_path, '_mibps.txt');
right_mibpfile = strcat(RNA_NAME, '_', right_tree_path, '_mibps.txt');
dlmwrite(left_mibpfile, new_mibps);
dlmwrite(right_mibpfile, new_mibps);
left_pfsname = strcat(RNA_NAME, '_', left_tree_path, '.pfs');
right_pfsname = strcat(RNA_NAME, '_', right_tree_path, '.pfs');
left_ctname = strcat(RNA_NAME, '_', left_tree_path, '.ct');
right_ctname = strcat(RNA_NAME, '_', right_tree_path, '.ct');


%make left child constraint file
if length(find(left_tree_path=='1')) > 0
   left_forced_pairs = new_mibps(left_tree_path == '1', :);
else
    left_forced_pairs = [];
end
if length(find(left_tree_path=='0')) > 0
   left_forced_nonpairs = new_mibps(left_tree_path == '0', :);
else
    left_forced_nonpairs = [];
end
left_filename = strcat(RNA_NAME, '_', left_tree_path, '.CON');
left_forced_pairs = [left_forced_pairs; FORCED_PAIRS];
left_forced_nonpairs = [left_forced_nonpairs; FORCED_NONPAIRS];
create_constraint_file([], [], left_forced_pairs, left_forced_nonpairs, left_filename);

%make right child constraint file
if length(find(right_tree_path=='1')) > 0
   right_forced_pairs = new_mibps(find(right_tree_path == '1'), :);
else
   right_forced_pairs = [];
end
if length(find(right_tree_path=='0')) > 0
   right_forced_nonpairs = new_mibps(find(right_tree_path == '0'), :);
else
    right_forced_nonpairs = [];
end
right_filename = strcat(RNA_NAME, '_', right_tree_path, '.CON');
right_forced_pairs = [right_forced_pairs; FORCED_PAIRS];
right_forced_nonpairs = [right_forced_nonpairs; FORCED_NONPAIRS];
create_constraint_file([], [], right_forced_pairs, right_forced_nonpairs, right_filename);

% create left child's pfs and ct files
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PARTITION_PROG ' -c ' left_filename ' ' seqfile ' ' left_pfsname ' >/dev/null']);

% create right child's pfs and ct files
system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PARTITION_PROG ' -c ' right_filename ' ' seqfile ' ' right_pfsname ' >/dev/null']);

end

