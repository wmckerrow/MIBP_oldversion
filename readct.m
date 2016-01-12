function [ samples ] = readct(ctfile, seq_length)
% Get column 5 from the ct file. Output is a NUM_SAMPLES x seq_length
% array samples(i,j) is pair of base j in sample i or 0 if it is
% unpaired.

program_constants;
%disp('in ct2dots. input ctfile = ');
%disp(ctfile);

samples = zeros(NUM_SAMPLES,seq_length);

fid = fopen(ctfile);

%disp('after calling fid = fopen(ctfile) we have fid = ');
%disp(fid);
for j = 1:NUM_SAMPLES
    fgetl(fid);
    for i = 1:seq_length
        line = fgetl(fid);
        C = regexp(line, ...
                   '(\d+)\s\w\s+\d+\s+\d+\s+(\d+)\s+\d+', ...
                   'tokens');
        if size(C)
            pos1 = str2num(cell2mat(C{1}(1)));
            pos2 = str2num(cell2mat(C{1}(2)));
            samples(j,pos1)=pos2;
        end
    end
end
fclose(fid);

end

