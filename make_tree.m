function [tree_strings,all_tree_strings] = make_tree_hibp(tree_strings,all_tree_strings, num_MIBPs_left, mi_cutoff, energy)
% This function makes a posterior-space tree with the greedy algorithm that
% splits the cluster that has the MIBP with the highest probability-mass-weighted average
% mutual information. It (and other functions it calls) creates all the
% files necessary to store information about the tree
% INPUTS: tree_strings, a cell array of strings representing the current
%           leaf clusters of the tree. I.e., tree_strings are the nodes
%           that could be split by picking out the next MIBP. When the
%           tree is first created, tree_strings is set to {''}, since the
%           tree has only one node.
%         num_MIBPs, an integer, the number of most informative base pairs
%           to pick out, a.k.a. the number of splits to perform.
%         energy, the total energy of the entire ensemble. Used to
%           calculate the probability mass of a cluster
% OUTPUT: tree_strings, a cell array of strings representing the leaf
%           clusters of the finished tree. e.g., tree_strings could equal
%           {'01', '00', '110', '111', '100', '101'}.

program_constants;

num_MIBPs_left

if num_MIBPs_left > 0
    
    %find next MIBP, called "mibp":
    
    max_MI = -inf; % initialize variables. max_MI is the maximum weighted mutual information of a cluster in tree_strings
    mibp = [];
    split_node = '';
    split_node_index = -1;
    for i = 1:length(tree_strings)
        tree_path = tree_strings{i};
        MI_file = strcat(RNA_NAME, '_', tree_path, '_MI.txt');
        if exist(MI_file,'file')
            MI = dlmread(MI_file); % get mutual information, MIBP of tree_path's cluster
        else
            disp('Calculating MI for ')
            tree_path
            pfsname = strcat(RNA_NAME, '_', tree_path, '.pfs');
            ctname = strcat(RNA_NAME, '_', tree_path, '.ct');
            system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ;' STOCHASTIC_PROG ' ' pfsname ' ' ctname ' >/dev/null']);
            MI=[0,0,0];
            [MI(1),MI(2:3)] = B_get_MI(tree_path);
            dlmwrite(MI_file, MI);
            MI(1)
        end
        mutual_info = MI(1);
        base_pair = MI(2:3);
        mutual_info = MI(1);
        base_pair = MI(2:3);
        % weight mutual_info by probability mass of cluster. we pick MIBP
        % with highest weighted_mi:
        weighted_mi = mutual_info * B_get_probmass2(tree_path, energy);
        if weighted_mi > max_MI % if weighted_mi > max_MI, set mibp to be base_pair and update variables
            mibp = base_pair;
            split_cluster = tree_path;
            split_node_index = i;
            max_MI = weighted_mi;
        end
    end
    max_MI
    if max_MI < mi_cutoff
        return
    end
    
    [tree_strings,all_tree_strings] = full_split(split_cluster,tree_strings,all_tree_strings,mibp);
        
    [tree_strings,all_tree_strings] = make_tree(tree_strings, all_tree_strings, num_MIBPs_left - 1,mi_cutoff, energy); % recursively call B_make_trees with one less mibp
end
end
