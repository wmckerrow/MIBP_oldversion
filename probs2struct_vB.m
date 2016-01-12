function s = probs2struct_vB(probs, n, Q, num_structs, included_bps, excluded_bps, conflicting_bps, tree_path, split_num)
% This function converts the following data into a matlab structure:
% probability (probs), sequence length (n),
% partition functions (Q, two values for with and without entropy constraints),
% num_structs (num_structs, two values), BPs forced present (included_bps),
% BPs forced absent (excluded_bps), and tree path (tree_path). The
% structures is then converted into a .json file by savejson.m.

s = struct('nucleotides', 1:n);
s.Q = Q;
s.num_structs = num_structs;
s.tree_path = tree_path;
s.split_num = num2str(split_num);
s.probs = {};
for i = 1:size(probs, 1)
    s.probs{end+1} = struct('nuc1', probs(i, 1), 'nuc2', probs(i, 2), 'p', probs(i, 3));
end

s.included = {};
for i = 1:size(included_bps, 1)
   s.included{end+1} = struct('nuc1', included_bps(i, 1), 'nuc2', included_bps(i,2), 'p', 2);
end
s.excluded = {};
for i = 1:size(excluded_bps, 1)
   s.excluded{end+1} = struct('nuc1', excluded_bps(i, 1), 'nuc2', excluded_bps(i,2), 'p', -1);
end
s.conflicting = {};
for i = 1:size(conflicting_bps, 1)
   s.conflicting{end+1} = struct('nuc1', conflicting_bps(i, 1), 'nuc2', conflicting_bps(i,2), 'p', -1);
end

end

