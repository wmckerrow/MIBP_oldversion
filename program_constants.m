pg_file = fopen('program_constants.txt','r');

line = fgetl(pg_file);

while ischar(line)
    eval(line);
    line = fgetl(pg_file);
end

fclose(pg_file);