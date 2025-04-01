# Comparing common strategies for ortholog selection used in phylogenomics

These codes are for constructing core BUSCOs on metazoan phylogeny.

Basically, all the scripts would work on the output of the OrthoFinder directories, by doing MAFFT-Linsi sequences alignment and IQ-TREE2 tree building, we use a customized script to identify the out-paralog OGs and in-paralog sequences, and filtering both for further analysis.


## S1-Data collection

12 chromosome level genome assemblies (10 metazoans and 2 outgroups), and these species are Amphimedon_queenslandica, Mnemiopsis_leydi, Spizellomyces_punctatus, Capitella_teleta, Ciona_intestinalis, Aplysia_californica, Acropora_millepora, Hoilungia_hongkongensis, Priapulus_caudatus, Tribolium_castaneum, Salpingoeca_dolichothecata, and Strongylocentrotus_purpuratus) to perform orthologs inference using BUSCO and OF with default parameters. 

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](https://myoctocat.com/assets/images/base-octocat.svg)
this screenshot is the picture of tree plus completeness

BUSCO: Metazoan 954 families
OrthoFinder: Orthogroup_Sequences 71249 OGs, Single_Copy_Orthologur_Sequences

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](https://myoctocat.com/assets/images/base-octocat.svg)
this is the screen shot from slides of original BUSCO and OrthoFinder output



## S2-Make comparisons

### S2.1-Copy the Orthogroup_Sequences directory to the working directory.

1. Using the bash script to group all OGs into different categories:

+ for the Orthogroup_Sequences dir, we have 71249 OGs, we then group them into 2:
- folder OF_1_3: 57459 OGs
  these are OGs with less than 3 sequences, cannot be used to build trees, so we have to take them out to another folder.
- folder OF_ALL: 13790 OGs
  

2. MAFFT alignment
   using the script below
3. IQtree tree building
   using the script below
  
### S2.2topology analysis (identify paralogs, i.e. in-paralogs and out-paralogs)
1. identify_paralogs_2.py
```
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 10:54:20 2024

@author: jesuslozanofernandez
Modified on July by Mingzhu Yang
"""


from ete3 import Tree
from collections import defaultdict
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python identify_paralogs.py <treefile>")
    sys.exit(1)

my_tree = sys.argv[1]
print(f"Input file: {my_tree}")

# 获取文件名和路径/Get file name and path
base_name = os.path.basename(my_tree)
file_name, file_ext = os.path.splitext(base_name)

# 加载树/ load trees
tree = Tree(my_tree, format=1)
print("Tree successfully loaded")

my_dict = defaultdict(list)

# 构建叶子节点字典/ Build dictionaries of leafs
for leaf in tree.get_leaf_names():
    my_dict[leaf.split('_')[0]].append(leaf)

inparalogs = []
outparalogs = []

# 识别paralogs/ identify paralogs
for k, v in my_dict.items():
    if len(v) > 1:
        paralogs = tuple(v)
        is_monophyletic = tree.check_monophyly(paralogs, "name", unrooted=True)
        if is_monophyletic[0]:
            inparalogs.append(str(paralogs) + ': In-Paralogs')
        else:
            outparalogs.append(str(paralogs) + ': Out-Paralogs')


inparalogs_file = f"{file_name}_inparalogs.txt"
outparalogs_file = f"{file_name}_outparalogs.txt"

with open(inparalogs_file, "w") as f:
    for line in inparalogs:
        f.write(line + "\n")

with open(outparalogs_file, "w") as f:
    for line in outparalogs:
        f.write(line + "\n")

print(f"In-Paralogs written to {inparalogs_file}")
print(f"Out-Paralogs written to {outparalogs_file}")
```
   
```
Input file: /user/work/qw23953/6_Compare_Software/8_Redo/OG0000000.fa.mafft.treefile
Tree successfully loaded
In-Paralogs written to OG0000000.fa.mafft_inparalogs.txt
Out-Paralogs written to OG0000000.fa.mafft_outparalogs.txt
Input file: /user/work/qw23953/6_Compare_Software/8_Redo/OG0000001.fa.mafft.treefile
Tree successfully loaded
In-Paralogs written to OG0000001.fa.mafft_inparalogs.txt
Out-Paralogs written to OG0000001.fa.mafft_outparalogs.txt
Input file: /user/work/qw23953/6_Compare_Software/8_Redo/OG0000002.fa.mafft.treefile
Tree successfully loaded
In-Paralogs written to OG0000002.fa.mafft_inparalogs.txt
Out-Paralogs written to OG0000002.fa.mafft_outparalogs.txt

```
2. Prepare the file of outparalog_family.txt
   It should looks like
```
OG0000000.fa
OG0000002.fa
OG0000003.fa
OG0000004.fa
OG0000005.fa
OG0000006.fa
OG0000007.fa
OG0000008.fa
OG0000011.fa
...
```
3. Copy all the file into a new folder called 3_all_OG_no_outparalogs/
   Move out-paralog files to target_folder="${folder}outparalogs_file/"
   using the script of move_outparalog_file.sh
```
#!/bin/bash

# 定义文件路径
folder="3_all_OG_no_outparalogs/"
move_list="outparalog_family.txt"
log_file="move_files.log"
target_folder="${folder}outparalogs_file/"

# 清空日志文件（如果已存在）
> "$log_file"

# 创建目标文件夹（如果不存在）
mkdir -p "$target_folder"

# 读取要移动的文件名列表
while IFS= read -r filename
do
    # 构造完整的文件路径
    filepath="${folder}${filename}"
    
    # 检查文件是否存在，如果存在则移动并记录日志
    if [ -f "$filepath" ]; then
        echo "Moving $filepath to $target_folder" | tee -a "$log_file"
        mv "$filepath" "$target_folder"
    else
        echo "File $filepath does not exist" | tee -a "$log_file"
    fi
done < "$move_list"

echo "File moving process complete. Check $log_file for details."
```


4. Further filtering the in-paralog sequences



   
6.   


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



BUSCO	OrthoFinder
All OGs (contains 4-more-multi)	vary	13264
All single OGs (contains 4-12 single-copy)	954	1659
Single-copy OGs (12 single-copy present)	211	365
