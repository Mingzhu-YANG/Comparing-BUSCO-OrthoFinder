import os

# Specify the path to the BUSCO and Orthofinder folders
busco_folder = '/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/3_BUSCO_2025_processed/final'
orthofinder_folder = '/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/copy_OF_4_N_list/of_4_n_8139/new_OG_1_N'

# Create output file
output_file_path = 'OF_1_n_results.txt'

# Open the output file and prepare to write the comparison results
with open(output_file_path, 'w') as output_file:
    output_file.write("busco_file_name\tSeq_ID\torthofinder_file_name\tSeq_ID\tlength of sites identitical\n")
    
    # Loop through each fasta file in the Orthofinder folder
    for orthofinder_file in os.listdir(orthofinder_folder):
        if orthofinder_file.endswith(".fa"):  # Make sure only process fasta file
            orthofinder_filepath = os.path.join(orthofinder_folder, orthofinder_file)
            
            # Extract the file name of the Orthofinder file as an identifier
            orthofinder_filename = orthofinder_file.split('.')[0]
            
            # Read the sequence information in the Orthofinder file
            orthofinder_sequences = {}
            with open(orthofinder_filepath, 'r') as orthofinder_handle:
                current_sequence = ""
                for line in orthofinder_handle:
                    if line.startswith(">"):
                        current_sequence = line.strip().split(" ")[0]  # Extract sequence ID
                        orthofinder_sequences[current_sequence] = ""
                    else:
                        orthofinder_sequences[current_sequence] += line.strip()
            
            # Loop through each fasta file in the BUSCO folder
            for busco_file in os.listdir(busco_folder):
                if busco_file.endswith(".faa"):  # Only fasta file
                    busco_filepath = os.path.join(busco_folder, busco_file)
                    
                    # Extract the file name of the Orthofinder file as an identifier
                    busco_filename = busco_file.split('.')[0]
                    
                    # Read the sequence information in the BUSCO file
                    busco_sequences = {}
                    with open(busco_filepath, 'r') as busco_handle:
                        current_sequence = ""
                        for line in busco_handle:
                            if line.startswith(">"):
                                current_sequence = line.strip().split(" ")[0]  # Extract sequence ID
                                busco_sequences[current_sequence] = ""
                            else:
                                busco_sequences[current_sequence] += line.strip()
                    
                    # Check if each sequence in the Orthofinder file exists in each sequence in the BUSCO file
                    for orthofinder_id, orthofinder_seq in orthofinder_sequences.items():
                        for busco_id, busco_seq in busco_sequences.items():
                            if orthofinder_seq in busco_seq:
                                # Write the result to a txt file
                                output_file.write(f"{busco_filename}\t{busco_id}\t{orthofinder_filename}\t{orthofinder_id}\t{len(orthofinder_seq)}\n")

print("Comparison results have been written to", output_file_path)
