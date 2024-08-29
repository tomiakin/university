% Specify the output file names
output_file_a = 'angles_a.txt';
output_file_b = 'angles_b.txt';

% Open the file for writing (angles_a.txt)
fid_a = fopen(output_file_a, 'w');

% Check if the file was opened successfully
if fid_a == -1
    error('Error opening file for writing');
end

% Open the file for writing (angles_b.txt)
fid_b = fopen(output_file_b, 'w');

% Check if the file was opened successfully
if fid_b == -1
    error('Error opening file for writing');
end

% Specify the CSV file name
% file_name = 'motor_angles_4.5mm_circle.csv';
file_name = 'motor_angles_102mm_triangle_finalFFFF.csv';
% file_name = 'motor_angles_102mm_triangle.csv';

% Read the entire data from the CSV file
data = csvread(file_name);

% Extract the 3rd and 4th columns
A = data(:, 3);
B = data(:, 4);

% Write the values of A to the file with commas after each value
fprintf(fid_a, 'a\n');
for i = 1:length(A)
    fprintf(fid_a, '%d,\n', A(i));
end

% Close the file for A
fclose(fid_a);

disp(['Data written to ' output_file_a]);

% Write the values of B to the file with commas after each value
fprintf(fid_b, 'b\n');
for i = 1:length(B)
    fprintf(fid_b, '%d,\n', B(i));
end

% Close the file for B
fclose(fid_b);

disp(['Data written to ' output_file_b]);
