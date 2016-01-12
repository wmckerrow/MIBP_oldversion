function bp_probs = read_bp_probs(filename)
% This function reads the probability file created by RNAstructure's
% ProbabilityPlot function (see
% http://rna.urmc.rochester.edu/Text/ProbabilityPlot.html, note that we use
% the -t flag). It converts this file into a matrix of base pair
% probabilities.
% INPUT: filename, a string, the name of the probability file that we're
%          reading
% OUTPUT: bp_probs, an m x 3 matrix. The first two entries of a row are two
%           nucleotides. The third entry of the row is the probability that
%           these two nucleotides form a base pair

bp_probs = [];
fid = fopen(filename, 'r');
ignore = fgetl(fid);
ignore = fgetl(fid); %throw away first two lines of file
line = fgetl(fid);

% this code is hacky. you might just want to trust that it works or rewrite
% it
while (line ~= -1)
   base1 = '';
   base2 = '';
   prob = '';
   on = 1;
   for i = 1:length(line)
       if length(strtrim(line(i))) == 0
           %if the character is a space, move on to next value
           on = on + 1;
       else
           %else keep adding digit to current value
           if (on == 1)
               base1 = strcat(base1, line(i));
           elseif (on == 2)
               base2 = strcat(base2, line(i));
           else
               prob = strcat(prob, line(i));
           end
       end
   end
   base1 = str2double(base1);
   base2 = str2double(base2);
   prob = str2double(prob);
   prob = 10^(-prob); %convert prob from -log10(probability) to probability
   bp_probs(end+1,:) = [base1, base2, prob];
   
   line = fgetl(fid);
end

fclose(fid);

end

