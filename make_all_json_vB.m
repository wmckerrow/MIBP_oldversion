function make_all_json_vB(strings)
% Creates one .json file per vertex in the tree. These files store the base
% pair probabilities, the MIBPs, the partition function, and number of
% structures
% INPUT: strings, a cell array of strings, where each string is a
%          "tree_path" of a vertex in the tree. e.g. '' is the root node, '0' is the
%           root node's left child, etc. (see B_master_script for details)

program_constants;

for i = 1:length(strings)
    tree_path = strings{i};
    disp('tree path:');
    disp(tree_path);
    probs_file = strcat(RNA_NAME, '_', tree_path, '_probs.txt');
    bp_probs = read_bp_probs(probs_file);
    if length(tree_path) > 0 % if vertex is not the root, then there are MIBPs
        mibp_file = strcat(RNA_NAME, '_', tree_path, '_mibps.txt');
        mibps = dlmread(mibp_file);
    else
        mibps = [];
    end
    qstruct_file = strcat(RNA_NAME, '_', tree_path, '_QandStruct.txt');
    qstruct = dlmread(qstruct_file);
    Q = qstruct(1, :);
    num_structs = qstruct(2, :);
    included = mibps(tree_path == '1', :);
    excluded = mibps(tree_path == '0', :);
    conflicting = [];
    split_num = num2str(floor(find(ismember(strings,strcat(tree_path,'0')), 1 )/2));
    s = probs2struct_vB(bp_probs, RNA_LENGTH, Q, num_structs, included, excluded, conflicting, tree_path, split_num);
    json_file = strcat(RNA_NAME, '_', tree_path, '.json');
    json = savejson('', s, json_file);
end

end

