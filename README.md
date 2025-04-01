# Comparing common strategies for ortholog selection used in phylogenomics

These codes are for constructing core BUSCOs on metazoan phylogeny


## S1-Data collection

12 chromosome level genome assemblies (10 metazoans and 2 outgroups), and these species are Amphimedon_queenslandica, Mnemiopsis_leydi, Spizellomyces_punctatus, Capitella_teleta, Ciona_intestinalis, Aplysia_californica, Acropora_millepora, Hoilungia_hongkongensis, Priapulus_caudatus, Tribolium_castaneum, Salpingoeca_dolichothecata, and Strongylocentrotus_purpuratus) to perform orthologs inference using BUSCO and OF with default parameters. 

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](https://myoctocat.com/assets/images/base-octocat.svg)
this screenshot is the picture of tree plus completeness

BUSCO: Metazoan 954 families
OrthoFinder: Orthogroup_Sequences 71249 OGs, Single_Copy_Orthologur_Sequences

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](https://myoctocat.com/assets/images/base-octocat.svg)
this is the screen shot from slides of original BUSCO and OrthoFinder output



## S2-Make comparisons

### Copy the Orthogroup_Sequences directory to the working directory.

1. Using the bash script to group all OGs into 4 categories:
- folder_1_3
- folder_4_12_single
- folder_4_12_multi
- folder_more_than_12

the script used would be 1_more_12.sh, 2_find_dup.sh

```
#!/bin/bash

# Source dir
source_directory="/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/multi_gene_families2/more_than_12"

# Target dir
target_directory="/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/multi_gene_families2/more_than_12/z_should_be_0"

# Create target dir if not exist
mkdir -p "$target_directory"

# Go through all fasta files
for fasta_file in "$source_directory"/*.fa; do
    # Get file names
    filename=$(basename "$fasta_file")
    # Get number of sequences in each fasta file
    sequence_count=$(grep -c '^>' "$fasta_file")
    
    # If sequences number more than 12, the move file to target dir
    if [ $sequence_count -gt 12 ]; then
        mv "$fasta_file" "$target_directory/$filename"
        echo "Moved $filename to $target_directory"
    fi
done
```
```
#!/bin/bash

# Source dir
input_directory="/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/multi_gene_families2/more_than_12"

# Target dir
output_directory="/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/multi_gene_families2/more_than_12/uniq"

# Create target dir if not exist
mkdir -p "$output_directory"

# Go through all fasta files
for fasta_file in "$input_directory"/*.fa; do
    # Get file names
    filename=$(basename "$fasta_file")

    # Check if duplicate ID exists
    # 检查是否存在重复的物种 ID
    duplicate=$(grep -o '^>[A-Za-z]\{4\}' "$fasta_file" | sort | uniq -d)

    # 如果不存在重复的物种 ID，则移动文件到目标目录
    if [ -z "$duplicate" ]; then
        mv "$fasta_file" "$output_directory/$filename"
        echo "Moved $filename to $output_directory"
    fi
done


```



main comparison script compare.py
```
import os

# Specify the path to the BUSCO and Orthofinder folders
busco_folder = '/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/BUSCO/busco_extracted'
orthofinder_folder = '/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/duplicate_four_to_12'

# Create output file
output_file_path = '/Users/qw23953/Desktop/PhD_project/Metazoa/compare_software/BUSCO_OF_duplicate_4_12.txt'

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
                if busco_file.endswith(".fasta"):  # Only fasta file
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

```


2. 






  ### 1_GET_SEQUENCES.py

 ```python


def group_by_gene_family(input_dir, output_dir):
    """
    Group sequences by gene family and save into separate .faa files.
    """
    os.makedirs(output_dir, exist_ok=True)
    gene_family_dict = defaultdict(list)

    for file_name in os.listdir(input_dir):
        if file_name.endswith(".faa"):
            file_path = os.path.join(input_dir, file_name)

            for record in SeqIO.parse(file_path, "fasta"):
                match = re.search(r'(\d+)at33208', record.id)
                if match:
                    family_id = match.group(1)
                    gene_family_dict[family_id].append(record)
                else:
                    print(f"Skipping sequence with unexpected format: {record.id}")

    for family_id, sequences in gene_family_dict.items():
        output_file = os.path.join(output_dir, f"{family_id}at33208.faa")
        with open(output_file, "w") as f:
            SeqIO.write(sequences, f, "fasta")

    print("FASTA files have been grouped by gene family and saved in the output directory.")

def main(base_dir, intermediate_dir, final_output_dir, keep_intermediates=False):
    """
    Main function to process .faa files: rename, merge, and group by gene family.
    """
    # Step 1: Rename sequences
    rename_sequences(base_dir)

    # Step 2: Merge sequences by species
    merge_sequences(base_dir, intermediate_dir)

    # Step 3: Group sequences by gene family
    group_by_gene_family(intermediate_dir, final_output_dir)

    # Clean up intermediate files if not needed
    if not keep_intermediates:
        for file in os.listdir(intermediate_dir):
            os.remove(os.path.join(intermediate_dir, file))
        print("Intermediate files have been removed.")

# Paths and settings
BASE_DIR = "/path/to/folder"
INTERMEDIATE_DIR = "/path/to/folder/intermediate"
FINAL_OUTPUT_DIR = "/path/to/folder/final"
KEEP_INTERMEDIATES = False  # Change to True if intermediate files should be kept

# Run the pipeline
main(BASE_DIR, INTERMEDIATE_DIR, FINAL_OUTPUT_DIR, KEEP_INTERMEDIATES)
 ```




BUSCO	OrthoFinder
All OGs (contains 4-more-multi)	vary	13264
All single OGs (contains 4-12 single-copy)	954	1659
Single-copy OGs (12 single-copy present)	211	365
