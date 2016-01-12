function cluster = which_cluster3(ctfile,all_tree_strings)
% This function finds the region that contains that native structure. If
% the output ends in '0' then it is not a retained cluster.
% INPUTS: A ct file for the native structure, the list of all node in the
% tree.

program_constants;

structure = read_one_ct(ctfile,RNA_LENGTH);
cluster = '';
inside = (sum(strcmp([cluster '1'],all_tree_strings))>0); % inside is TRUE when cluster is not a leaf

while inside
    mibp_file = strcat (RNA_NAME, '_', cluster, '1_mibps.txt');
    mibp = dlmread(mibp_file);
    mibp = mibp(end,:);
    if structure(mibp(1)) == mibp(2) || structure(mibp(2)) == mibp(1)
        cluster = strcat(cluster,'1');
    else
        cluster = strcat(cluster,'0');
    end
    inside = (sum(strcmp([cluster '1'],all_tree_strings))>0);
end

end
