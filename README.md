# Comparing common strategies for ortholog selection used in phylogenomics

These codes are for constructing core BUSCOs on metazoan phylogeny.

Basically, all the scripts would work on the output of the OrthoFinder directories, by doing MAFFT-Linsi sequences alignment and IQ-TREE2 tree building, we can get the topologies of different orthogroups, then we use a customized script to identify the out-paralog OGs and in-paralog sequences(here we only keep one copy to make them single copy, notice we only keep the one that busco have kept), and filtering both for further analysis.


## S1-Data collection

12 chromosome level genome assemblies (10 metazoans and 2 outgroups, and these species are Amphimedon_queenslandica, Mnemiopsis_leydi, Spizellomyces_punctatus, Capitella_teleta, Ciona_intestinalis, Aplysia_californica, Acropora_millepora, Hoilungia_hongkongensis, Priapulus_caudatus, Tribolium_castaneum, Salpingoeca_dolichothecata, and Strongylocentrotus_purpuratus) to perform orthologs inference using BUSCO and OrthoFinder with default parameters. 

![12 chromosome level genome assemblies](/final_combined_plot_reversed_fixed2.png)
this screenshot is showing the taxa we selected and their completeness

BUSCO: Metazoan 954 families

OrthoFinder: Orthogroup_Sequences 71249 OGs, Single_Copy_Orthologur_Sequences

![Screenshot of a comment on a GitHub issue showing an image, added in the Markdown, of an Octocat smiling and raising a tentacle.](/busco_of_dir_output.png)
this is the screenshot of original BUSCO and OrthoFinder output

------------------------------------------------------------------------------------------------------------------------------------------------------------------


## S2-Make comparisons
### S2.0- Get all the busco families using get_sequences.py (remember to create a copy before change)

Go to the output path of BUSCO, it should looks like:

```
qw23953@it103273:~/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025$ ls
Acropora_millepora.fasta.out        Mnemiopsis_leydi.fasta.out
Amphimedon_queenslandica.fasta.out  Priapulus_caudatus.fasta.out
Aplysia_californica.fasta.out       Salpingoeca_dolichothecata.fasta.out
busco_downloads                     Spizellomyces_punctatus.fasta.out
Capitella_teleta.fasta.out          Strongylocentrotus_purpuratus.fasta.out
Ciona_intestinalis.fasta.out        Tribolium_castaneum.fasta.out
Hoilungia_hongkongensis.fasta.out
```

***get_sequences.py***

```
import os
from Bio import SeqIO
import re
from collections import defaultdict

def rename_sequences(base_dir):
   """
   Rename sequences in .faa files under .out folders by appending folder and file names.
   """
   for folder in os.listdir(base_dir):
       if folder.endswith(".out"):
           folder_name = folder.replace(".fasta.out", "")
           faa_dir = os.path.join(base_dir, folder, "run_metazoa_odb12", "busco_sequences", "single_copy_busco_sequences")

           if not os.path.exists(faa_dir):
               print(f"Directory not found: {faa_dir}")
               continue

           for faa_file in os.listdir(faa_dir):
               if faa_file.endswith(".faa"):
                   faa_path = os.path.join(faa_dir, faa_file)
                   faa_name = faa_file.replace(".faa", "")

                   new_records = []
                   for record in SeqIO.parse(faa_path, "fasta"):
                       new_id = f"{folder_name}_{faa_name}"
                       record.id = new_id
                       record.description = ""
                       new_records.append(record)

                   SeqIO.write(new_records, faa_path, "fasta")
                   print(f"Processed file: {faa_path}")

   print("All faa sequences are renamed!")

def merge_sequences(base_dir, output_dir):
   """
   Merge renamed .faa files into a single file for each species (folder).
   """
   os.makedirs(output_dir, exist_ok=True)

   for folder in os.listdir(base_dir):
       if folder.endswith(".out"):
           source_dir = os.path.join(base_dir, folder, "run_metazoa_odb12/busco_sequences/single_copy_busco_sequences")
           if not os.path.exists(source_dir):
               print(f"Warning: {source_dir} does not exist, skipping.")
               continue

           output_filename = f"{folder.replace('.fasta.out', '')}.faa"
           output_file = os.path.join(output_dir, output_filename)

           with open(output_file, "w") as outfile:
               for faa_file in os.listdir(source_dir):
                   if faa_file.endswith(".faa"):
                       faa_path = os.path.join(source_dir, faa_file)
                       for record in SeqIO.parse(faa_path, "fasta"):
                           SeqIO.write(record, outfile, "fasta")

           print(f"Created {output_file}.")

   print("All faa files are merged!")

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
BASE_DIR = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025_copy"
INTERMEDIATE_DIR = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025_copy/intermediate"
FINAL_OUTPUT_DIR = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/1_BUSCO_2025_copy/final"
KEEP_INTERMEDIATES = True  # Change to True if intermediate files should be kept

# Run the pipeline
main(BASE_DIR, INTERMEDIATE_DIR, FINAL_OUTPUT_DIR, KEEP_INTERMEDIATES)
```
By using the script above, you will get two new folders in your working directory, they are :

   -/path/to/intermediate
   
   -/path/to/final

You will need the data from /final folder, these are all the busco single-copy gene families, contains all the species sequences belong to that single-copy genes.
later on you would use this dataset to run comparisons.
   
------------------------------------------------------------------------------------------------------------------------------------------------------------------



### S2.1-Copy the Orthogroup_Sequences directory to the working directory (now we work on OrthoFinder output).

1. Using the bash script to group all OGs into different categories, for the Orthogroup_Sequences dir, we have 71249 OGs, we then group them into 2:
   
   - folder OF_1_3: 57459 OGs, using script **get_groups_1_3_seqs.sh**
     these are OGs with less than 3 sequences, cannot be used to build trees, so we have to take them out to folder OF_1_3_seqs.
   - folder OF_4_N: 13790 OGs, we separate them using script **find_single_copy.sh**
       - we further select the OGs with 4-12 sequences without duplication, there are 1659 OGs
       - so the rest 12131 OGs with duplications
       - we can also use **get_strict_single_copy.sh** to get single-copy orthogroups with all taxa present

**get_groups_1_3_seqs.sh**
```
#!/bin/bash

# Define directory
SOURCE_DIR="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451"   
TARGET_DIR="$SOURCE_DIR/OF_1_3_seqs"

# Create target dir (if not exist)
mkdir -p "$TARGET_DIR"

# Loop through all the fa file
for file in "$SOURCE_DIR"/*.fa; do
    # Caculate the number of seqs in each file
    count=$(grep -c "^>" "$file")

    # Check number of seqs, whether they are in between the defined number
    if [ "$count" -ge 1 ] && [ "$count" -le 3 ]; then
        # Move file to target folder
        mv "$file" "$TARGET_DIR/"
        # echo "移动文件: $file"
    fi
done

echo "Done！"
```
**find_single_copy.sh**
```
#!/bin/bash

input_directory="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/OF_4_N"

output_directory="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/OF_4_N/OF_4_12_single"

mkdir -p "$output_directory"

# Loop through all the fasta file
for fasta_file in "$input_directory"/*.fa; do
    # Get the file name
    filename=$(basename "$fasta_file")
    # Check if there are duplicated species ID
    duplicate=$(grep -o '^>[A-Za-z]\{4\}' "$fasta_file" | sort | uniq -d)
    # If no duplicated ID, then move file to target folder
    if [ -z "$duplicate" ]; then
        mv "$fasta_file" "$output_directory/$filename"
        echo "Moved $filename to $output_directory"
    fi
done
```
**get_strict_single_copy.sh**

```
#!/bin/bash

# Define directory
SOURCE_DIR="/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/OF_4_N/OF_4_12_single"   
TARGET_DIR="$SOURCE_DIR/OF_12_seqs"

# Create target dir (if not exist)
mkdir -p "$TARGET_DIR"

# Loop through all the fa file
for file in "$SOURCE_DIR"/*.fa; do
    # Caculate the number of seqs in each file
    count=$(grep -c "^>" "$file")

    # Check number of seqs, whether they are in between the defined number
    if [ "$count" -ge 12 ]; then
        # Move file to target folder
        mv "$file" "$TARGET_DIR/"
        echo "Moved file: $file"
    fi
done
```




2. MAFFT alignment

```
nohup bash -c 'for i in *.fa; do mafft --localpair --maxiterate 1000 --preservecase "$i" > "$i.mafft"; done' > mafft.log 2>&1 &
```

       
3. IQtree tree building
   
    - 196 OGs are failed from MAFFT or IQtree, so copy all the successful .fa and .treefile to a new folder called 2_all_OG_multi_13594 (or 11935)

```
for i in *.mafft
do

iqtree2 -s $i -m MFP -mset LG+F+G -madd LG+F+G+C60,LG+F+G+C50,LG+C40+F+G,LG+C30+F+G,LG+C20+F+G,LG+C10+F+G,C10,C20,C30,C40,C50,C60 --score-diff all -mwopt  -bb 1000

done
```
 Or if there are too many files or too big files, you need to configure the parallel tasks to avoid disrupt the system
```
find . -maxdepth 1 -iname "*mafft" > job_array.txt
```
```
nohup bash -c '
    nt=4
    max_parallel_jobs=25
    cat job_array.txt | parallel -j $max_parallel_jobs --no-notice "
        echo \"Running {}\"
        iqtree2 -s {} \
            -m MFP \
            -mset LG+F+G \
            -madd LG+F+G+C60,LG+F+G+C50,LG+C40+F+G,LG+C30+F+G,LG+C20+F+G,LG+C10+F+G,C10,C20,C30,C40,C50,C60 \
            --score-diff all \
            -mwopt \
            -nt ${nt} \
            -bb 1000 > {}.log 2>&1
    "
' > parallel.log 2>&1 &
```

----------------------------------------------------------------------------------------------------------------------------------------------------------------

 
### S2.2-topology analysis (identify paralogs, i.e. in-paralogs and out-paralogs)
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

# Get file name and path
base_name = os.path.basename(my_tree)
file_name, file_ext = os.path.splitext(base_name)

# Load trees
tree = Tree(my_tree, format=1)
print("Tree successfully loaded")

my_dict = defaultdict(list)

# Build dictionaries of leafs
for leaf in tree.get_leaf_names():
    my_dict[leaf.split('_')[0]].append(leaf)

inparalogs = []
outparalogs = []

# Identify paralogs
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

# Check if files are empty and delete if so
for file in [inparalogs_file, outparalogs_file]:
    if os.path.getsize(file) == 0:
        os.remove(file)
        print(f"Deleted empty file: {file}")
```
 the output should like below:  
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
OG0000000_tree_outparalogs.txt
OG0000001_tree_outparalogs.txt
OG0000002_tree_outparalogs.txt
OG0000003_tree_outparalogs.txt
OG0000004_tree_outparalogs.txt
...
```
3. Copy all the file into a new folder called 3_all_OG_no_outparalogs/
   
   Move out-paralog files to target_folder="${folder}outparalogs_file/"
   
   using the script of move_outparalog_file.sh
```
import os
import shutil

# Base directory containing all the files
base_dir = "/home/qw23953/mingzhu/4_BUSCO_OrthoFinder_2025/4_Orthogroup_Sequences_63451/copy_OF_4_N_list/of_4_n_8139"
# Target directory to move the matched files into
output_dir = os.path.join(base_dir, "outparalog_family")
# Text file that contains the list of *_outparalogs.txt filenames
list_file = os.path.join(base_dir, "outparalog_family.txt")

# Create the output directory if it doesn't already exist
os.makedirs(output_dir, exist_ok=True)

# Read the list file line by line
with open(list_file, "r") as f:
    for line in f:
        line = line.strip()
        if line.endswith("_outparalogs.txt"):
            # Remove the "_outparalogs.txt" suffix to get the prefix
            prefix = line.replace("_outparalogs.txt", "")
            # Construct the corresponding tree file name
            tree_file = f"{prefix}.txt"
            src_path = os.path.join(base_dir, tree_file)
            dst_path = os.path.join(output_dir, tree_file)
            
            # Move the file if it exists
            if os.path.exists(src_path):
                shutil.move(src_path, dst_path)
                print(f"✅ Moved: {tree_file}")
            else:
                print(f"❌ File not found: {tree_file}")        
```
and move all the *_outparalogs.txt file and *inparalogs.txt file by edit tree_file = f"{prefix}.txt" line of codes, the suffix part.

4. Further filtering the in-paralog sequences
   
   Here we going to filter the in-paralog sequences, we keep only once copy thus make it single copy, we keep the sequence that busco have kept.

   To do this, we first cat all the inparalog file (e.g. OG0000000.fa.mafft_inparalogs.txt) into one file called inparalogs_after_out.txt

   it looks like
   
```
('HHON_TR10104_c1_g1_i1_m_21166', 'HHON_TR22145_c0_g1_i1_m_43391',  'HHON_TR2361_c1_g1_i2_m_2943', 'HHON_TR28545_c1_g1_i2_m_60573', 'HHON_TR28545_c1_g2_i1_m_60574'): In-Paralogs
('SDOL_gnl_est_initial1030.p2', 'SDOL_gnl_est_initial10959.p1', 'SDOL_gnl_est_initial8225.p2', 'SDOL_gnl_est_initial4818.p1', 'SDOL_gnl_est_initial26086.p2'): In-Paralogs
('HHON_TR10068_c4_g1_i2_m_20869', 'HHON_TR10068_c5_g1_i1_m_20870', 'HHON_TR22316_c3_g2_i2_m_43441', 'HHON_TR29456_c0_g1_i1_m_63671', 'HHON_TR146_c1_g1_i1_m_17'): In-Paralogs
('MLEY_Mnemiopsis_leydi_ML000128a', 'MLEY_Mnemiopsis_leydi_ML003255a', 'MLEY_Mnemiopsis_leydi_ML327420a', 'MLEY_Mnemiopsis_leydi_ML02332a'): In-Paralogs
('STRO_XP_011676996.2', 'STRO_XP_780641.3'): In-Paralogs
```
- 1. First we count the in-paralogs by using count_inparalogs.py
```
from collections import Counter

# Read file
with open("inparalogs_after_out.txt", "r") as file:
    data = file.readlines()

# Initialize the list
sequences = []

# Loop through each line to get the seq ID
for line in data:
    # Get the part in single quotation
    start = line.find("(")
    end = line.find(")")
    
    if start != -1 and end != -1:
        # Extract the sequences and remove extra symbol
        seqs = line[start+1:end].replace("'", "").split(", ")
        # Add to the seq list
        sequences.extend(seqs)

# Count the number of occurence of each seq ID
sequence_count = Counter(sequences)

# Print out seq ID and times
for seq, count in sequence_count.items():
    print(f"{seq}: {count}")

# Export all ID into file
with open("inparalogs_4_n_index.txt", "w") as output_file:
    for seq in sequences:
        output_file.write(f"{seq}\n")

# Pring number of sequences
print(f"Total number of sequences: {len(sequences)}")
```
 Then the output should be a file called  inparalogs_4_n_index.txt and it looks like:
```
HHON_TR10104_c1_g1_i1_m_21166
HHON_TR22145_c0_g1_i1_m_43391
HHON_TR3035_c1_g3_i10_m_6505
HHON_TR25907_c2_g1_i1_m_51961
HHON_TR14192_c0_g2_i1_m_27580
HHON_TR28827_c0_g1_i1_m_61152
HHON_TR23333_c1_g1_i1_m_46165
```

sometimes need further process by using
```
awk '{
    split($0, a, "_");
    out = a[3];
    for (i = 4; i <= length(a); i++) {
        out = out "_" a[i];
    }
    print out;
}' inparalogs_4_n_index.txt > inparalogs_4_n_index_extracted.txt

```


- 2. Extract all the sequences above using extract_inparalogs_seq.py
```
from Bio import SeqIO
import os

# Path
index_file = "inparalogs_4_n_index_extracted.txt"
fasta_dir = "/user/work/qw23953/6_Compare_Software/8_Redo/4_all_OG_no_inparalogs/"
output_file = "inparalogs_4_n.fasta"

# Read the sequence list
with open(index_file, 'r') as f:
    target_sequences = set(line.strip() for line in f)

# Initialize the list
sequences_found = []
sequences_not_found = []

# Loop through all .fa file and extract seqs
for fasta_file in os.listdir(fasta_dir):
    if fasta_file.endswith(".fa"):
        file_path = os.path.join(fasta_dir, fasta_file)
        # Parse the FASTA file
        for record in SeqIO.parse(file_path, "fasta"):
            # If the sequence name is in the target sequence name list
            if record.id in target_sequences:
                sequences_found.append(record)
                target_sequences.remove(record.id)  # Removes the found sequence

# Save the extracted sequence to a new fasta file
with open(output_file, 'w') as output_handle:
    SeqIO.write(sequences_found, output_handle, 'fasta')

# Print the sequence name that was not found
sequences_not_found = list(target_sequences)
if sequences_not_found:
    print("Below sequence was not found:")
    for seq in sequences_not_found:
        print(seq)

# Print the total number of extracted sequences
print(f"Total number of extracted sequences: {len(sequences_found)}")
```

  Then we got an output of inparalogs_4_n.fasta, next step nwe need to map them with BUSCO sequences and keep them single copy, here we use 
  find_inparalogs_busco.py
  
```
from Bio import SeqIO
import os

# PATH
input_fasta = "inparalogs_4_n.fasta"
busco_dir = "/user/work/qw23953/6_Compare_Software/5_Sep_redo/busco_12_taxa"
output_file = "hits_match.txt"
to_delete_file = "inparalogs_to_delete.txt"

# Read 1736 sequences in inparalogs_4_n.fasta
query_sequences = {}
for record in SeqIO.parse(input_fasta, "fasta"):
    query_sequences[str(record.seq)] = record.id

# Initializes the matching record
matches = []
found_sequences = set()

# Iterate over all fasta files in the busco_12_taxa directory
for fasta_file in os.listdir(busco_dir):
    if fasta_file.endswith(".fasta"):
        file_path = os.path.join(busco_dir, fasta_file)
        for record in SeqIO.parse(file_path, "fasta"):
            seq_str = str(record.seq)
            if seq_str in query_sequences:
                # Record the matching result: busco file name and matching sequence name
                matches.append(f"{fasta_file}: {record.id}")
                found_sequences.add(seq_str)

# Generate the hits_match.txt file
with open(output_file, 'w') as f:
    for match in matches:
        f.write(match + "\n")

# Save the sequence of hits not found to inparalogs_to_delete.txt and count
not_found_count = 0
with open(to_delete_file, 'w') as f:
    for seq_str, seq_id in query_sequences.items():
        if seq_str not in found_sequences:
            f.write(f">{seq_id}\n{seq_str}\n")
            not_found_count += 1

# Count the number of sequences that find a match
found_count = len(found_sequences)
total_queries = len(query_sequences)

# Output statistics
print(f"In {total_queries} sequences, found {found_count} matchings ")
print(f"Matching results were save to {output_file} ")
print(f"Unmatching sequences were saved to{to_delete_file} ")
print(f"Number of unmatching sequences {not_found_count} ")
```
Then we got out oput of hits_match.txt and inparalogs_to_delete.txt

- 3. Delete all the inparalogs sequences above using delete_inparalogs.py
```
from Bio import SeqIO
import os

# Path
delete_fasta = "inparalogs_to_delete.txt"
target_dir = "/user/work/qw23953/6_Compare_Software/8_Redo/4_all_OG_no_inparalogs"

# Read the sequence to be deleted
sequences_to_delete = set()
for record in SeqIO.parse(delete_fasta, "fasta"):
    sequences_to_delete.add(str(record.seq))

# Count the number of deleted sequences
deleted_sequences_count = 0

# Loop through all .fa files in the directory
for fa_file in os.listdir(target_dir):
    if fa_file.endswith(".fa"):
        file_path = os.path.join(target_dir, fa_file)
        temp_file_path = file_path + ".temp"

        # Read the original file and writes to the temporary file, excluding sequences to be deleted
        with open(file_path, 'r') as infile, open(temp_file_path, 'w') as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                if str(record.seq) not in sequences_to_delete:
                    SeqIO.write(record, outfile, 'fasta')
                else:
                    deleted_sequences_count += 1

        # Replace original file
        os.replace(temp_file_path, file_path)

# Output
print(f"Number of deleted sequences is: {deleted_sequences_count}")
```
   - 4. Then, a lot of OrthoGroups become null file as they have no sequences remained, we move them to null_file folder using a script move_null_file.sh
   
```
#!/bin/bash

source_folder="/user/work/qw23953/6_Compare_Software/8_Redo/4_all_OG_no_inparalogs"
target_folder="${source_folder}/null_file"


mkdir -p "$target_folder"

# Find empty files and move to the destination folder
find "$source_folder" -type f -empty -exec mv {} "$target_folder" \;

echo "All empty files have been moved to $target_folder."
```


   
5. To now, for 13594 OGs, 3972 are out-paralog OGs, after removing in-paralog sequences, we lost 3993 OGs, so 5629 OGs are clean to use after filtering, they are in a folder called OF_4_N.
   
   + For now, we get the most relaxed level of single-copy Orthologues Groups.
   + Before we filtering, we also have OF_4_12 using BUSCO definition of single copy orthologs, the number of OF_4_12 is 1659 OGs,
   + We also have OrthoFinder definition of 365 Single_Copy_Orthologue_Sequences.

  **Using these diferent level of restrictions we make conparison with BUSCO output, here the BUSCO output was reformatted according to gene families.**
   

  main comparison script compare.py
```
import os

# Specify the path to the BUSCO and Orthofinder folders
busco_folder = '/user/work/qw23953/6_Compare_Software/4_Aug_redo/busco_12_taxa'
orthofinder_folder = '/user/work/qw23953/6_Compare_Software/8_Redo/5_all_filtered'

# Create output file
output_file_path = '/user/work/qw23953/6_Compare_Software/8_Redo/OF_all_results.txt'

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


