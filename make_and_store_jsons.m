program_constants
load([RNA_NAME '.mat'])
make_all_Q_struct_files_vB(all_tree_strings)
make_all_json_vB(all_tree_strings)
system(['mv -f ',RNA_NAME,'*.json jsons/']);

split_name = strsplit(RNA_NAME,'/');
system(['sed s/RNA_NAME_HERE/',split_name{end},'/ vis_template.html > temp1.html']);
vis_tree_string = strcat('\''',all_tree_strings{1},'\''');
for i = 2:length(all_tree_strings)
    vis_tree_string = strcat(vis_tree_string,',','\''',all_tree_strings{i},'\''');
end
system(['sed s/TREE_STRING_HERE/',vis_tree_string,'/ temp1.html > temp2.html']);
tree_depth = 0;
for i = 1:length(all_tree_strings)
    if tree_depth < length(all_tree_strings{i})
        tree_depth = length(all_tree_strings{i});
    end
end
system(['sed s/TREE_DEPTH_HERE/',num2str(tree_depth),'/ temp2.html > visualization.html']);