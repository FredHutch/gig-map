#!/usr/bin/env python3

import os

input_folder = 'inputs'
output_file = '${sample_name}.num_reads.txt'
total_sum = 0

# Iterate through each file in the input folder
for filename in os.listdir(input_folder):
    file_path = os.path.join(input_folder, filename)
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            try:
                number = int(file.readline().strip())
                total_sum += number
            except ValueError:
                print(f"Skipping file {filename} as it does not contain a valid integer.")

# Write the total sum to the output file
with open(output_file, 'w') as file:
    file.write(str(total_sum))
