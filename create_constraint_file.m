function create_constraint_file(forced_paired, forced_unpaired, forced_pairs, forced_nonpairs, file_name)
% This function creates a constraint file in the format specified by
% RNAstructure. (See
% http://rna.urmc.rochester.edu/Text/File_Formats.html#Constraint)
% This file specifies which nucleotides/pairs of nucleotides must or cannot
% be paired.
% INPUTS: forced_paired, an array, nucleotides that must be in a base pair.
%         forced_unpaired, an array, nucleotides that cannot be in a base pair.
%         forced_pairs, an n x 2 matrix, pairs of nucleotides that must
%           form a base pair.
%         forced_nonpairs, an n x 2 matrix, pairs of nucleotides that
%           cannot form a base pair.
%         file_name, a string, the name of the constraint file. Must end in
%           ".C"

fid = fopen(file_name, 'w');
fprintf(fid, 'DS:\n');
for i = 1:length(forced_paired)
    nuc = num2str(forced_paired(i));
    wrt = strtrim(strcat(nuc, '\n'));
    fprintf(fid, wrt);
end
fprintf(fid, '-1\nSS:\n');
for i = 1:length(forced_unpaired)
   nuc = num2str(forced_unpaired(i));
   wrt = strtrim(strcat(nuc, '\n'));
   fprintf(fid, wrt);
end
fprintf(fid, '-1\nMod:\n-1\nPairs:\n');
for i = 1:size(forced_pairs, 1)
   nuc1 = strtrim(num2str(forced_pairs(i, 1)));
   nuc2 = strtrim(num2str(forced_pairs(i, 2)));
   wrt = horzcat(nuc1, ' ', nuc2, '\n');
   fprintf(fid, wrt);
end
fprintf(fid, '-1 -1\nFMN:\n-1\nForbids:\n');
for i = 1:size(forced_nonpairs, 1)
   nuc1 = strtrim(num2str(forced_nonpairs(i, 1)));
   nuc2 = strtrim(num2str(forced_nonpairs(i, 2)));
   wrt = horzcat(nuc1, ' ', nuc2, '\n');
   fprintf(fid, wrt);
end
fprintf(fid, '-1 -1');
fclose(fid);

end

