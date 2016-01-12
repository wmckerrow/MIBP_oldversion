function make_all_probs_files_vB(strings)
% This function converts all the vertices' pfs files into probability (of
% base pairs) files
% INPUT: strings, a cell array of strings, where each string is a
%          "tree_path" of a vertex in the tree. e.g. '' is the root node, '0' is the
%           root node's left child, etc. (see B_master_script for details)


program_constants;

for i = 1:length(strings)
    disp(strings{i});
    tree_path = strings{i};
    pfsfile = strcat(RNA_NAME, '_', tree_path, '.pfs');
    probs_file = strcat(RNA_NAME, '_', tree_path, '_probs.txt');
    system(['export LD_LIBRARY_PATH=' GCC '; export DATAPATH=' DATAPATH ' ; ' PROB_PROG ' -t ' pfsfile ' ' probs_file ' >/dev/null']);
end

end

