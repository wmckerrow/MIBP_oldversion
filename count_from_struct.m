function counts = count_from_struct(samples)
% Count how many times each base pair appears in a sample.
% INPUT: A 1000 x RNA_LENGTH array describing which pairs apear in which
% samples
% OUTPUT: An RNA_LENGTH x RNA_LENGTH array with base pair counts.
counts = zeros(size(samples,2));
for i = 1:size(samples,1) % Loop through samples
    for j = 1:size(samples,2) % Loop through sequence
        pairedto = samples(i,j);
        if j < pairedto
            % Add to the array if there is a base pair.
            counts(j,pairedto) = counts(j,pairedto)+1;
        end
    end
end