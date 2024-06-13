import csv
import sys
import os

# Function to filter and save a TSV file
def generate_filtered_tsv(input_tsv, sample_name):
    filtered_tsv_file = f'{sample_name}.filtered.tsv'
    total_mapped_reads = 0

    # Check if the output file already exists
    file_exists = os.path.exists(filtered_tsv_file)

    # Open input TSV file
    with open(input_tsv, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        # Initialize a list to store filtered rows
        filtered_rows = []

        # Define the header for the output TSV file
        header = ["mapping_ref", "ref_length", "mapped_reads", "percentage_of_mapped_reads"]

        # Iterate through rows in the input TSV file
        for row in reader:
            # Check if the row has enough elements
            if len(row) >= 3:
                try:
                    # Extract the value from column 2 and add it to the total
                    value_to_check = int(row[2])
                    total_mapped_reads += value_to_check

                    # Append the row to the filtered rows list if the value is not 0
                    if value_to_check != 0:
                        filtered_rows.append(row[:3] + row[4:])
                except (IndexError, ValueError):
                    # Skip invalid rows and print a message
                    print(f"Skipped invalid row: {row}")
            else:
                print(f"Row does not have enough elements: {row}")

        # Calculate percentages for each row in the filtered rows list
        for row in filtered_rows:
            mapped_reads = int(row[2])
            percentage = (mapped_reads / total_mapped_reads) * 100
            row.append(f'{percentage:.2f}%')

    # Write the filtered data to the output TSV file
    with open(filtered_tsv_file, 'a' if file_exists else 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        if not file_exists:
            writer.writerow([f"Sample Name: {sample_name}"])
            writer.writerow(header)
        writer.writerows(filtered_rows)

    return filtered_tsv_file

# Function to generate refs TSV file
def generate_refs_tsv(filtered_tsv, sample_name, read_length, mean_depth_threshold):
    refs_tsv_file = f'{sample_name}.refs.tsv'
    with open(filtered_tsv, 'r') as infile, open(refs_tsv_file, 'w', newline='') as outfile:
        tsv_reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        next(tsv_reader)  # Skip the header
        for row in tsv_reader:
            try:
                if (int(row[2]) * read_length) / int(row[1]) >= mean_depth_threshold:
                    writer.writerow([row[0]])
            except (IndexError, ValueError):
                print(f"Skipped invalid row: {row}")
    return refs_tsv_file

# Entry point of the script
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 formater.py <input_tsv_file> <sample_name> <read_length> <mean_depth_threshold>")
        sys.exit(1)

    input_tsv = sys.argv[1]
    sample_name = sys.argv[2]
    read_length = int(sys.argv[3])
    mean_depth_threshold = int(sys.argv[4])

    # Check if input files exist
    for file in [input_tsv]:
        if not os.path.exists(file):
            print(f"Error: File not found - {file}")
            sys.exit(1)

    filtered_tsv = generate_filtered_tsv(input_tsv, sample_name)
    refs_tsv = generate_refs_tsv(filtered_tsv, sample_name, read_length, mean_depth_threshold)
    print(f"Generated filtered TSV file: {filtered_tsv}")
    print(f"Generated refs TSV file: {refs_tsv}")
