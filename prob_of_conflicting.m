for i = 1:length(all_tree_strings)
	cluster = all_tree_strings{i};
	if length(cluster) > 0 && cluster(end) == '0'
		temp = dlmread([RNA_NAME,'_',cluster,'_mibps.txt']);
		bp = temp(end,:);
		[conflicting_bp,prob] = get_conflicting_of_pairs(bp,cluster);
		disp(['conf_prob ',cluster,':',num2str(prob)])
		bps = [bp;conflicting_bp];
		for i = 1:4
			[conflicting_bp,newprob] = get_conflicting_of_pairs(bps,cluster);
			prob = prob+newprob;
			bps = [bps;conflicting_bp];
		end
		disp(['5conf_prob ',cluster,':',num2str(prob)])
	end
end